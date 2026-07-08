"""W3 Component D (data-prep half) — foam spherized-URDF parsing + the FLANGE anchor
mapping + self-collision ranges for `grid_collision`.

These are the pure-Python helpers that turn a foam-spherized URDF into a sphere batch
descriptor compatible with `_multitarget.build_target_batch` (W1b/W2a). They are
deliberately DECOUPLED from the codegen/cli/registry wiring (the `grid_collision`
namespace emitter + `--collision` flag co-land later with the W1b.3/W2a.3 registration),
so the ⚠FLANGE mapping — the #1 silent-wrong-frame hazard — can be unit-tested standalone
against real GRiD robots WITHOUT the (heavy, optional) foam toolchain.

KEY DIFFERENCE from HJCD's `foam_spheres.py`: GRiD binds each sphere to ITS OWN GRiD frame
(the `s_Xworld` joint slot), resolved through GRiD's URDFParser — NOT foam's actuated-joint
count and NOT HJCD's URDF-ordinal table. A sphere is orientation-free, so a sphere on a
WELDED (fixed-joint) frame reduces to its nearest MOVABLE ancestor joint with a constant
PRE-COMPOSED offset (`offset_parent = (T_parent<-welded @ [offset_local,1])[:3]`) — no
fixed-anchor `s_Xworld` slot required. `T_parent<-welded` is GRiD's OWN post-`remove_fixed_
joints` transform (`Fixed_Joint.get_transformation_matrix_hom()`), keeping GRiD's baked FK
the single source of truth.
"""
import xml.etree.ElementTree as ET

import numpy as np


# --------------------------------------------------------------------------- foam parse
def parse_spherized_urdf(path):
    """foam output: per link, one `<collision><geometry><sphere radius/></geometry>
    <origin xyz/></collision>` per sphere (center in the LINK frame; rpy irrelevant for a
    sphere). Returns `{link_name: [(x, y, z, radius), ...]}` in URDF document order."""
    root = ET.parse(path).getroot()
    out = {}
    for link in root.findall("link"):
        name = link.get("name")
        spheres = []
        for col in link.findall("collision"):
            sph = col.find("geometry/sphere")
            if sph is None:
                continue
            r = float(sph.get("radius"))
            origin = col.find("origin")
            xyz = (origin.get("xyz") if origin is not None else "0 0 0").split()
            x, y, z = (float(v) for v in xyz)
            spheres.append((x, y, z, r))
        if spheres:
            out[name] = spheres
    return out


def urdf_joint_tree(path):
    """Parse the URDF's native joint tree (retained in foam's spherized output). Returns
    `{child_link_name: joint_name}` — the bridge from a foam sphere's LINK to the GRiD
    Fixed_Joint (matched by joint name) when that link is welded and has no movable joint."""
    root = ET.parse(path).getroot()
    child_to_joint = {}
    for joint in root.findall("joint"):
        child = joint.find("child")
        if child is not None:
            child_to_joint[child.get("link")] = joint.get("name")
    return child_to_joint


# --------------------------------------------------------------------------- FLANGE mapping
def sphere_anchor_frames(robot, link_names, urdf_path):
    """⚠FLANGE crux. Map each foam link name -> `(anchor_jid, T_compose)` where `anchor_jid`
    is the GRiD movable-joint slot in `s_Xworld` and `T_compose` (4x4) carries the sphere's
    LOCAL offset into that anchor's frame (identity for a movable link; the fixed transform
    for a welded frame). A sphere at local `r` extracts at world `X_world[anchor] @ T_compose
    @ [r, 1]`.

    Resolution per link:
      * MOVABLE link  -> the joint whose child IS this link (`get_joints_by_child_name`);
        `T_compose = I`.
      * WELDED frame  -> the URDF joint feeding this link (by child-link name) matches a GRiD
        `Fixed_Joint` by NAME; anchor = that fixed joint's re-parented movable joint
        (`get_parent()` -> movable joint id), `T_compose = Fixed_Joint.get_transformation_
        matrix_hom()` (GRiD's post-removal parent-movable -> welded transform).
      * ROOT/WORLD frame (fed by NO joint at all, or a fixed joint welded to world parent -1)
        -> `(-1, I)` SKIP sentinel; caller drops these (pedestal-bolted) spheres.

    A mismatch here silently checks collision on the WRONG frame — validate against the
    UR10e/panda/iiwa base-0-monotone-down-chain assertion (see test_collision_flange_mapping)."""
    child_to_joint = urdf_joint_tree(urdf_path)
    frames = {}
    for link_name in link_names:
        movable = robot.get_joints_by_child_name(link_name)
        if movable:
            frames[link_name] = (int(movable[0].get_id()), np.eye(4))
            continue
        # No movable joint has this child. Bridge foam link -> URDF joint name -> GRiD Fixed_Joint.
        jname = child_to_joint.get(link_name)
        if jname is None:
            # No joint feeds this link -> it is the robot root/base frame. Skip (pedestal).
            frames[link_name] = (-1, np.eye(4))
            continue
        fj = robot.get_fixed_joint_by_name(jname)
        if fj is not None:
            parent_joint = robot.get_joint_by_name(fj.get_parent())
            if parent_joint is None:
                # fixed joint welded straight to the world/base (parent -1) -> skip sentinel
                frames[link_name] = (-1, np.eye(4))
                continue
            T = np.asarray(fj.get_transformation_matrix_hom(), dtype=np.float64).reshape(4, 4)
            frames[link_name] = (int(parent_joint.get_id()), T)
            continue
        raise ValueError(
            "FLANGE: cannot map foam link '%s' to a GRiD frame (no movable joint with this "
            "child, and no fixed joint named '%s'). A silent wrong-frame anchor would result."
            % (link_name, jname))
    return frames


def _parent_jid_map(robot):
    """`{jid: parent_jid}` over movable joints (parent movable joint = the joint whose child
    link is this joint's parent link; -1 at the root). Used for self-collision adjacency."""
    child_link_to_jid = {j.get_child(): j.get_id() for j in robot.get_joints_ordered_by_id()}
    parent = {}
    for j in robot.get_joints_ordered_by_id():
        parent[j.get_id()] = child_link_to_jid.get(j.get_parent(), -1)
    return parent


def _anchors_adjacent(parent_jid, a, b):
    """Two sphere anchors are 'adjacent' (skip the pair: they always touch) when they share a
    frame or are directly parent/child in the movable-joint tree."""
    return a == b or parent_jid.get(a, -1) == b or parent_jid.get(b, -1) == a


def build_self_cc_ranges(robot, anchor):
    """Self-collision pairs as `{sphere_i, start_j, end_j}` rows: sphere i is checked against
    spheres [start_j..end_j] (j > i). Adjacent-frame pairs (same/parent/child link) are
    skipped (always in contact). Same-anchor spheres are contiguous in `anchor` (built link by
    link), so the non-adjacent partners of i compress into maximal contiguous [j0,j1] runs."""
    parent_jid = _parent_jid_map(robot)
    n = len(anchor)
    ranges = []
    for i in range(n):
        j = i + 1
        while j < n:
            if _anchors_adjacent(parent_jid, anchor[i], anchor[j]):
                j += 1
                continue
            j0 = j
            while j < n and not _anchors_adjacent(parent_jid, anchor[i], anchor[j]):
                j += 1
            ranges.append((i, j0, j - 1))
    return ranges


def build_sphere_tiers(robot, foam_outputs):
    """`foam_outputs = {tier_name: spherized_urdf_path}` (e.g. broad + fine, two foam runs).
    Returns `{tier: {"n", "anchor"[N], "offset"[3N], "radius"[N], "self_cc_ranges"[R][3]}}`.
    Offsets are PRE-COMPOSED into their anchor frame (welded frames folded onto the movable
    parent). Index-0 (base/pedestal) spheres are skipped, matching HJCD."""
    tiers = {}
    for tier_name, path in foam_outputs.items():
        parsed = parse_spherized_urdf(path)
        frames = sphere_anchor_frames(robot, list(parsed.keys()), path)
        anchor, offset, radius = [], [], []
        for link_name, spheres in parsed.items():
            anchor_jid, T = frames[link_name]
            if anchor_jid < 0:
                continue  # skip root/world (pedestal-bolted) spheres
            for (x, y, z, r) in spheres:
                p = T @ np.array([x, y, z, 1.0])
                anchor.append(int(anchor_jid))
                offset.extend([float(p[0]), float(p[1]), float(p[2])])
                radius.append(float(r))
        tiers[tier_name] = {
            "n": len(anchor), "anchor": anchor, "offset": offset, "radius": radius,
            "self_cc_ranges": build_self_cc_ranges(robot, anchor),
        }
    return tiers
