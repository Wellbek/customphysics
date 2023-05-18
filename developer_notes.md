Developer notes
===============

This document describes how this module was constructed. Students should not need to
know these details. It is intended to facilitate maintenance (e.g. Godot version upgrades).

The purpose of the `customphysics` module is to create an alternative to the Bullet module
in which the actual physics is left out, but collision detection and overall integration
with Godot is left intact.

The main steps taken to produce this module are summarized below.

1. Copy all .cpp and .h files from the bullet module into the `godot_module_files` subdirectory,
   except the `register_types.*` files.
2. Search-and-replace in `godot_module_files`: `Bullet` -> `Custom`.
   This requires special care. In general, we do not want to do this in
   any include statements, we only want to rename *types introduced by the
   Godot bullet integration*. However, it seems that Bullet itself only
   uses the "bt" prefix, and only uses "Bullet" in the file/directory
   names for includes, so a direct search-and-replace *almost* works. Using a tool
   that gives a good manual overview is necessary, however.
   This is necessary in order to ensure that we can compile `customphysics` in the same
   code base as the original `bullet` module, as otherwise we'll run into name clash problems,
   and moreover there are some two-way connections between types that prevent us
   from reusing more of the `bullet` infrastructure.
3. Rename the following funtions in the files `godot_module_files/bullet_types_converter.*`:
     - B_TO_G -> B_TO_G_CUSTOM
     - G_TO_B -> G_TO_B_CUSTOM
     - INVERT_B_TO_G -> INVERT_B_TO_G_CUSTOM
     - INVERT_G_TO_B -> INVERT_G_TO_B_CUSTOM
     - UNSCALE_BT_BASIS -> UNSCALE_BT_BASIS_CUSTOM
   This can be accomplished with a search-and-replace as before.
4. To prevent name clashes, the following functions must also be renamed:
     - `calculateGodotCombinedRestitution` -> `calculateGodotCombinedRestitution_custom`
     - `calculateGodotCombinedFriction` -> `calculateGodotCombinedFriction_custom`
     - `equal` -> `equal_custom` (local function in collision_object_bullet.cpp)
5. In `space_bullet.cpp`, the function `create_empty_world` must be modified to create
   instances of `CustomDynamicsWorld` instead of Bullet's dynamic world types.
6. The following types must be renamed accordingly:
    - btRayShape
    - GodotAllContactResultCallback
    - GodotAllConvexResultCallback
    - GodotClosestConvexResultCallback
    - GodotClosestRayResultCallback
    - GodotCollisionConfiguration
    - GodotCollisionDispatcher
    - godotContactAddedCallback
    - GodotContactPairContactResultCallback
    - GodotDeepPenetrationContactResultCallback
    - GodotFilterCallback
    - GodotKinClosestConvexResultCallback
    - GodotRayWorldAlgorithm
    - GodotRestInfoContactResultCallback
    - GodotSoftCollisionConfiguration
