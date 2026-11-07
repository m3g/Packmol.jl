# rigid_body

The molecules are rigid: their variables are the center of mass and three Euler angles.

| File | Contents |
|---|---|
| `rigid_body.jl` | Rotation matrices (`eulermat`), moving and rotating molecules (`move!`, `rotate!`, `random_move!`), the alignment of a molecule to observed positions (`kabsch_rotation`, `euler_angles`), and the reference coordinates of the structure types (`set_reference_coordinates!`: free molecules centered, with the longest axis along z). |
| `chain_rule.jl` | The gradient relative to the variables of the molecules, from the gradient relative to the atom positions (`chain_rule!`). |
| `rotation_metric.jl` | The scaling of the rotation variables by the inertia of each molecule (`RotationScaling`): the optimizer works on scaled angles, a unit step of which moves the atoms by 1 Å RMS. |
