# # # # # # # # # # # # # # # # # # # # # # # # #
#      _____                   __  __           #
#     |_   _| __ _   _ ___ ___|  \/  | ___      #
#       | || '__| | | / __/ __| |\/| |/ _ \     #
#       | || |  | |_| \__ \__ \ |  | |  __/     #
#       |_||_|   \__,_|___/___/_|  |_|\___|     #
#                                               #
# A Maple Library for Truss Elements Structures #
# # # # # # # # # # # # # # # # # # # # # # # # #

# Authors: Matteo Larcher and Davide Stocco
# Date:    21/12/2022

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

with(plots): # FIXME: better not to include and to specify plot in the ffunction
# call (e.g. "plots:-plot3d(...)")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# TODO: create a 'TrussMe' module
#TrussMe := module()
#export FRAME, GROUND, BEAM, ROD, FORCE, MOMENT, QLOAD, SUPPORT, JOINT, MATERIAL,
#       STRUCTURE, Show, Rotate, Translate, Project, DefineMaterial, DefineStructure,
#       DefineBeam, DefineRod, DefineJoint, DefineSupport, DefineForce,
#       DefineMoment, DefineQload, Solve;
#end module:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   _____
#  |_   _|   _ _ __   ___  ___
#    | || | | | '_ \ / _ \/ __|
#    | || |_| | |_) |  __/\__ \
#    |_| \__, | .__/ \___||___/
#        |___/|_|

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

`type/FRAME`     := {Matrix}:
`type/GROUND`    := {table}:
`type/BEAM`      := {table}:
`type/ROD`       := {table}:
`type/FORCE`     := {table}:
`type/MOMENT`    := {table}:
`type/QLOAD`     := {table}: # FIXME: Change name for 'QLOAD' in something like 'QFORCE', 'DIST_FORCE' or 'DISTRUBUTED_FORCE'
#`type/QTORQUE`     := {table}: # TODO: Consider introducing distributed torque
`type/SUPPORT`   := {table}:
#`type/COMPLIANTSUPPORT`   := {table}: # TODO: Consider introducing compliant joint
`type/JOINT`     := {table}:
`type/MATERIAL`  := {table}:
`type/STRUCTURE` := {table}:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   ___       _                        _
#  |_ _|_ __ | |_ ___ _ __ _ __   __ _| |
#   | || '_ \| __/ _ \ '__| '_ \ / _` | |
#   | || | | | ||  __/ |  | | | | (_| | |
#  |___|_| |_|\__\___|_|  |_| |_|\__,_|_|
#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ground := <<1, 0, 0, 0>|
           <0, 1, 0, 0>|
           <0, 0, 1, 0>|
           <0, 0, 0, 1>>:

_gravity := [0, 0, 0]:

EARTH := table({
  type             = 'GROUND',
  name             = "ground",
  length           = 0,
  frame            = ground,
  admissible_loads = [1, 1, 1, 1, 1, 1]
  }):

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   ____            _            _
#  |  _ \ _ __ ___ | |_ ___  ___| |_
#  | |_) | '__/ _ \| __/ _ \/ __| __|
#  |  __/| | | (_) | ||  __/ (__| |_
#  |_|   |_|  \___/ \__\___|\___|\__|
#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# FIXME: capitalize or not? (e.g. `type` or `TYPE`), do we need to add/remove someone more?
protect(
  'type',
  'ground',
  'name',
  'length',
  'frame',
  'admissible_loads',
  'elastic_modulus',
  'poisson_modulus',
  'shear_modulus',
  'density',
  'components',
  'coordinate',
  'force',
  'moment',
  'forces',
  'moments',
  'beam',
  'qload',
  'dof_constr',
  'cs_inertias',
  'target',
  'cs_area',
  'internal_actions',
  'material',
  'X',
  'Y',
  'Z',
  'TWOD'
  );

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   _____                 _   _
#  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___
#  | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
#  |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
#  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Show := proc(
  tab::{table}, # Table to be shown
  $)

  description "Show the content of a table";

  print(tab = tab[type](op(op(tab))));
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

GetNames := proc(
  objs::{ # List/set of structural elements
    list({'MATERIAL', 'BEAM', 'ROD', 'SUPPORT', 'JOINT'}),
    set( {'MATERIAL', 'BEAM', 'ROD', 'SUPPORT', 'JOINT'})
  }, $)

  description "Get names of a list/set of structural elements";

  return [seq(objs[i][name], i = 1..nops(objs))];
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Rotate := proc(
  axis,  # Rotation axis
  angle, # Rotation angle
  $)

  description "Rotation around an axis 'X','Y' or 'Z' of a given angle";

  if (axis = 'X') then
    return <<1, 0,           0,          0>|
            <0, cos(angle),  sin(angle), 0>|
            <0, -sin(angle), cos(angle), 0>|
            <0, 0,           0,          1>>;
  elif (axis = 'Y') then
    return <<cos(angle), 0, -sin(angle), 0>|
            <0,          1, 0,           0>|
            <sin(angle), 0, cos(angle),  0>|
            <0,          0, 0,           1>>;
  elif (axis = 'Z') then
    return <<cos(angle),  sin(angle), 0, 0>|
            <-sin(angle), cos(angle), 0, 0>|
            <0,           0,          1, 0>|
            <0,           0,          0, 1>>;
  else
    error "wrong axis detected";
  end
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Translate := proc(
  x, # X-axis translation component
  y, # Y-axis translation component
  z, # Z-axis translation component
  $)

  description "Translation on x-, y- and z-axis";

  return <<1, 0, 0, 0>|
          <0, 1, 0, 0>|
          <0, 0, 1, 0>|
          <x, y, z, 1>>;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Project := proc(
  x,       # Vector/point to be projected
  RF_from, # Reference frame from which the vector/point is expressed
  RF_to,   # Reference frame to which the vector/point will be expressed
  $)

  description "Project vector[x,y,z,0]/point[x,y,z,1] x from RF_from to RF_to";

  LinearAlgebra[MatrixInverse](RF_to).RF_from.<x[1], x[2], x[3], x[4]>;
  return simplify([seq(%[i],i = 1..4)]);
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   __  __       _            _       _
#  |  \/  | __ _| |_ ___ _ __(_) __ _| |
#  | |\/| |/ _` | __/ _ \ '__| |/ _` | |
#  | |  | | (_| | ||  __/ |  | | (_| | |
#  |_|  |_|\__,_|\__\___|_|  |_|\__,_|_|
#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

DefineMaterial := proc({
  id::{string} := "steel",      # Name of the material
  E            := 210.0E+09,    # Elastic modulus (Pa)
  nu           := 0.3,          # Poisson modulus (-)
  G            := E/(2*(1+nu)), # Shear modulus (Pa)
  rho          := 7.4E+03       # Density (kg/m^3)
  }, $)

  description "Define a 'MATERIAL' object with inputs: name of the material, "
    "elastic modulus E (default = 210.0E9 Pa), Poisson modulus nu (default = "
    "0.3), shear modulus G (default = E/(2*(1+nu))), density rho (default = "
    "7.4E3 kg/m^3)";

  return table({
    name            = id,
    type            = 'MATERIAL',
    elastic_modulus = E,
    poisson_modulus = nu,
    shear_modulus   = G,
    density         = rho
    });
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   _____
#  |  ___|__  _ __ ___ ___
#  | |_ / _ \| '__/ __/ _ \
#  |  _| (_) | | | (_|  __/
#  |_|  \___/|_|  \___\___|
#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeForce := proc(
  comps::{list},                  # Force components
  ell,                            # Application point (axial coordinate)
  obj::{'BEAM', 'ROD', 'GROUND'}, # Target object
  RF::{'FRAME'} := ground,        # Reference frame in which the force is defined
  $)

  description "Define a 'FORCE' object with inputs: force components, force "
    "application axial coordinate [0,L], target object, optional reference "
    "frame in which the force is defined (default = ground)";

  local proj_comps;

  if (obj[type] = 'BEAM') or (obj[type] = 'ROD') then
    if (ell < 0) or (ell > obj[length]) then
      error "force application point must be in [0,L] range";
    end if;
  end if;

  proj_comps := Project([op(comps), 0], RF, obj[frame])[1..3];
  if (obj[type] = 'ROD') then
    if (proj_comps[2] <> 0) or (proj_comps[3] <> 0) then
      error "only axial forces are accepted in 'ROD' objects";
    end if;
  elif (obj[type] = 'SUPPORT') or (obj[type] = 'JOINT') then
    if (ell <> 0) then
      error "only null axial coordinate is accepted for 'SUPPORT' and 'JOINT' "
        "objects";
    end if;
  end if;

  return table({
    type       = 'FORCE',
    components = proj_comps,
    coordinate = ell,
    target     = obj[name]
    });
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeMoment := proc(
  comps::{list},                  # Moment components
  ell,                            # Application point (axial coordinate)
  obj::{'BEAM', 'ROD', 'GROUND'}, # Target object
  RF::{'FRAME'} := ground,        # Reference frame in which the moment is defined
  $)

  description "Define a 'MOMENT' object with inputs: moment components, "
    "moment application axial coordinate [0,L], target object, optional "
    "reference frame in which the moment is defined (default = ground)";

  local proj_comps;

  if (obj[type] = 'BEAM') or (obj[type] = 'ROD') then
    if (ell < 0) or (ell > obj[length]) then
      error "moment application point must be in [0,L] range";
    end if;
  end if;

  proj_comps := Project([op(comps),0],RF,obj[frame])[1..3];
  if (obj[type] = 'ROD') then
    "moment cannot be applied to 'ROD' objects";
    # FIXME: error "moment cannot be applied to 'ROD' objects";
  elif (obj[type] = 'SUPPORT') or (obj[type] = 'JOINT') then
    if (ell <> 0) then
      error "only null axial coordinate is accepted for 'SUPPORT' and 'JOINT' "
        "objects";
    end if;
  end if;

  return table({
    type       = 'MOMENT',
    components = proj_comps,
    coordinate = ell,
    target     = obj[name]
    });
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeQLoad := proc(
  comps::{list},                # Distributed load components
  obj::{'BEAM', 'ROD'},         # Target object
  ell_min       := 0,           # Initial application point (axial coordinate)
  ell_max       := obj[length], # Final application point (axial coordinate)
  RF::{'FRAME'} := ground,      # Reference frame in which the object is defined
  $)

  description "Define a 'QLOAD' object with inputs: load components, target "
    "object, initial, and final application points (axial coordinates), "
    "optional reference frame in which the load components are defined "
    "(default = ground)";

  local proj_comps;

  proj_comps := Project([op(comps), 0], RF, obj[frame])[1..3];
  if (obj[type] = 'ROD') then
    if (proj_comps[2] <> 0) or (proj_comps[3] <> 0) then
      error "only axial loads are accepted in 'ROD' objects"
    end if;
  end if;

  return table({
    type       = 'QLOAD',
    components = [
      x -> piecewise((x >= ell_min) and (x <= ell_max), proj_comps[1], 0),
      x -> piecewise((x >= ell_min) and (x <= ell_max), proj_comps[2], 0),
      x -> piecewise((x >= ell_min) and (x <= ell_max), proj_comps[3], 0)
    ],
    target     = obj[name]
    });
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# FIXME: consider deleting the support object and merge it with the joint object
MakeSupport := proc(
  id::{string},                # Support name
  dof_constrained::{list},     # Constrained degree of freedom
  objs::list({'BEAM', 'ROD'}), # Target objects
  ells::{list},                # Support locations
  RF::{'FRAME'} := ground,     # Reference frame in which the support is defined
  $)

  description "Define a 'SUPPORT' object with inputs: support name, constrained degree of freedom, target objects, list of support locations, optional reference frame in which the support is defined (default = ground)";

  local S, J_tmp, i, j, sr;
  global EARTH;

  for i from 1 to nops(objs) do
    if (objs[i][type] = 'ROD') and (ells[i] <> 0) and (ells[i] <> objs[i][length]) then
      error "'SUPPORT' objects can only be applied at extremes of 'ROD' objects"
    end if;
    if (objs[i][type] = 'ROD') and (dof_constrained[4..6] <> [0, 0, 0]) then
      error "'ROD' objects supports can only have translational constraints"
    end if;
  end do;

  S := table({
    type                     = 'SUPPORT',
    dof_constr               = dof_constrained,
    coords                   = [0, op(ells)],
    name                     = id,
    frame                    = RF,
    targets                  = [EARTH[name], op(GetNames(objs))],
    variables                = [],
    forces                   = [],
    moments                  = [],
    constraint_loads         = [],
    constraint_displacements = [],
    support_reactions        = [] # Expressed in ground reference frame
    });

  # Build the temporary joint
  J_tmp := MakeJoint(id, dof_constrained, [EARTH, op(objs)], S[coords], RF);

  S[variables]                := J_tmp[variables];
  S[forces]                   := J_tmp[forces];
  S[moments]                  := J_tmp[moments];
  S[constraint_loads]         := J_tmp[constraint_loads];
  S[constraint_displacements] := J_tmp[constraint_displacements];

  # Retrieve support force reactions
  sr := [FX, FY, FZ, MX, MY, MZ];
  for i from 1 to nops(S[forces]) do
    if (S[forces][i][target] = "ground") then
      for j from 1 to 3 do
        if (S[forces][i][components][j] <> 0) then
          S[support_reactions] := [
            op(S[support_reactions]),
            sr[j] = -S[forces][i][components][j]
            ];
        end if;
      end do;
      break;
    end if;
  end do;

  # Retrieve support moments reactions
  for i from 1 to nops(S[moments]) do
    if (S[moments][i][target] = "ground") then
      for j from 1 to 3 do
        if (S[moments][i][components][j] <> 0) then
          S[support_reactions] := [
            op(S[support_reactions]),
            sr[j+3] = -S[moments][i][components][j]
            ];
        end if;
      end do;
      break;
    end if;
  end do;

  return op(S);
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanSupport := proc(
  obj::{'SUPPORT'}, # Support to be cleaned
  $)

  description "Clean 'SUPPORT' object internal variables";

  # TODO: check if this is necessary
  #obj[variables]                := [];
  #obj[forces]                   := [];
  #obj[moments]                  := [];
  #obj[constraint_loads]         := [];
  #obj[constraint_displacements] := [];
  #obj[support_reactions]        := [];
  return obj;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeJoint := proc(
  id::{string},                # Joint name
  dof_constrained::{list},     # Constrained degree of freedom
  objs::list({'BEAM', 'ROD'}), # Target objects
  X::{list},                   # Joint locations
  RF::{'FRAME'} := ground,     # Reference frame in which the joint is defined
  $)

  description "Make a 'JOINT' object with inputs: joint name, constrained "
    "degrees of freedom, target objects, joint locations, and optional "
    "reference frame in which the joint is defined (default = ground)";

  local J, i, jf_comp, jm_comp, jf_comp_obj, jm_comp_obj, jm_indets, jf_indets,
    constraint;
  global EARTH;

  for i from 1 to nops(objs) do
    if (objs[i][type] = 'ROD') and (X[i] <> 0) and (X[i] <> objs[i][length]) then
      error "'JOINT' objects can only be applied at extremes of 'ROD' objects";
    end if;
    if (objs[i][type] = 'ROD') and (dof_constrained[4..6] <> [0, 0, 0]) then
      error "'ROD' objects supports can only have translational constraints";
    end if;
  end do;


    J := table({
      type = 'JOINT',
      dof_constr               = dof_constrained,
      coords                   = X,
      name                     = id,
      frame                    = RF,
      targets                  = GetNames(objs),
      variables                = [],
      forces                   = [],
      moments                  = [],
      constraint_loads         = [],
      constraint_displacements = []
      });

  # Add all the bodies forces
  for i from 1 to nops(objs) do
    # Create force compatible with the joint constrained dof
    jf_comp := convert(<
      JFx_||(J[name])||_||(objs[i][name]),
      JFy_||(J[name])||_||(objs[i][name]),
      JFz_||(J[name])||_||(objs[i][name])
      > *~ <op(dof_constrained[1..3])>,
      list);
    # Project the components into object frame and extract admissible loads
    jf_comp_obj := convert(
      Project([op(jf_comp), 0], RF, objs[i][frame])[1..3]
      .~ <op(objs[i][admissible_loads][1..3])>,
      list);
    # Use the non admissible loads to build the loads constraint
    constraint := convert(
      Project([op(jf_comp), 0], RF, objs[i][frame])[1..3]
      .~ <op((-1*objs[i][admissible_loads][1..3]) +~ 1)>,
      list);
    constraint := remove(x -> x = 0, constraint);
    J[constraint_loads] := [
      op(J[constraint_loads]),
      op(constraint)
      ];
    # Extract the survived components
    jf_indets := indets(jf_comp);
    # Check if there are reactions
    if (jf_comp_obj <> [0, 0, 0]) then
      # Create the reaction force between joint and obj
      JF_||(id)||_||(objs[i][name]) := MakeForce(jf_comp_obj, X[i], objs[i], objs[i][frame]);
      JF_||(objs[i][name])||_||(id) := MakeForce(-jf_comp_obj, 0, J, objs[i][frame]);
      # Update the output joint
      J[variables] := [
        op(J[variables]),
        op(jf_indets)
        ];
      J[forces] := [
        op(J[forces]),
        JF_||(id)||_||(objs[i][name]),
        JF_||(objs[i][name])||_||(id)
        ];
    end if;
  end do;

  # Add all the bodies moments
  for i from 1 to nops(objs) do
    # Create moment compatible with joint constrained dof
    jm_comp := convert(<
      JMx_||(J[name])||_||(objs[i][name]),
      JMy_||(J[name])||_||(objs[i][name]),
      JMz_||(J[name])||_||(objs[i][name])
      > *~ <op(dof_constrained[4..6])>,
      list);
    # Project the components into object frame and extract the admissible loads
    jm_comp_obj := convert(
      Project([op(jm_comp), 0], RF, objs[i][frame])[1..3]
      .~ <op(objs[i][admissible_loads][4..6])>,
      list);
    # Use the non admissible loads to build the loads constraint
    constraint := convert(
      Project([op(jm_comp), 0], RF, objs[i][frame])[1..3]
      .~ <op((-1*objs[i][admissible_loads][4..6]) +~ 1)>,
      list);
    constraint := remove(x -> x = 0, constraint);
    J[constraint_loads] := [
      op(J[constraint_loads]),
      op(constraint)
      ];
    # Extract the survived components
    jm_indets := indets(jm_comp);
    # Check if there are reactions
    if (jf_comp_obj <> [0, 0, 0]) then
      # Create the reaction force between joint and obj
      JM_||(id)||_||(objs[i][name]) := MakeMoment(jm_comp_obj, X[i], objs[i], objs[i][frame]);
      JM_||(objs[i][name])||_||(id) := MakeMoment(-jm_comp_obj, 0, J, objs[i][frame]);
      # Update the output joint
      J[variables] := [
        op(J[variables]),
        op(jm_indets)
        ];
      J[moments] := [
        op(J[moments]),
        JM_||(id)||_||(objs[i][name]),
        JM_||(objs[i][name])||_||(id)
        ];
    end if;
  end do;

  return op(J);
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanJoint := proc(
  obj::{'JOINT'}, # Object to be cleaned
  $)

  description "Clean 'JOINT' object internal variables";

  # TODO: check if this is necessary
  #obj[internal_actions] := [];
  return obj;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeRod := proc(
  id::{string}, # Object name
  RF,           # Reference frame
  ell,          # Length (m)
  {
    A                 := 0,   # Cross-section area (m^2)
    mat::{'MATERIAL'} := NULL # Material
  },
  $)

  description "Create a 'ROD' object with inputs: object name, reference "
  "frame, length, and optional cross-section area and material";

  return table({
    type             = 'ROD',
    name             = id,
    length           = ell,
    cs_area          = A,
    material         = mat,
    frame            = RF,
    admissible_loads = [1, 0, 0, 0, 0, 0],
    internal_actions = []
    });
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanRod := proc(
  obj::{'ROD'}, # Object to be cleaned
  $)

  description "Clean 'ROD' object internal variables";

  # TODO: check if this is necessary
  #obj[internal_actions] := [];
  return obj;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeBeam := proc(
  id::{string}, # Object name
  RF,           # Reference frame
  ell,          # Length (m)
  {
    A                 := 0,   # Cross-section area (m^2)
    I_xx              := 0,   # Cross-section x-axis inertia (m^4)
    I_yy              := 0,   # Cross-section y-axis inertia (m^4)
    I_zz              := 0,   # Cross-section z-axis inertia (m^4)
    mat::{'MATERIAL'} := NULL # Material
  },
  $)

  description "Create a 'BEAM' object with inputs: object name, reference "
    "frame, length, and optional cross-section area, inertias on x-, y- and "
    "z-axis and material";

  return table({
    type             = 'BEAM',
    name             = id,
    length           = ell,
    cs_area          = A,
    cs_inertias      = [I_xx, I_yy, I_zz],
    material         = mat,
    frame            = RF,
    admissible_loads = [1, 1, 1, 1, 1, 1],
    internal_actions = []
    });
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanBeam := proc(
  obj::{'BEAM'}, # Object to be cleaned
  $)

  description "Clean 'BEAM' object internal variables";

  # TODO: check if this is necessary
  #obj[internal_actions] := [];
  return obj;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeStructure := proc(
  objs::{ # Structure objects (BEAM,ROD,SUPPORT,JOINT)
    list({'BEAM', 'ROD', 'SUPPORT', 'JOINT'}),
    set( {'BEAM', 'ROD', 'SUPPORT', 'JOINT'})
  },
  ext::{ # External forces, moments or distributed loads
    list({'FORCE', 'MOMENT', 'QLOAD'}),
    set( {'FORCE', 'MOMENT', 'QLOAD'})
  } := [],
  {
    hyper_vars::{list ,set} := [],
      # Hyperstatic variables
    hyper_disp::{list, set} := [seq(0, 1..nops(hyper_vars))],
      # Hyperstatic displacements
    dim::{string} := "3D"
      # Structure dimension ("2D" or "3D")
  },
  $)

  description "Create a 'STRUCTURE' object with inputs: structure objects, "
    "external forces, moments or distributed loads, hyperstatic variables and "
    "displacements and structure dimension (""2D"" or ""3D"")";

  local num_dof, i, names, candidate_hyp_vars;

  # Check for duplicate names
  names := [];
  for i from 1 to nops(objs) do
    if member(objs[i][name], names) then
      error "duplicate names found on structure objects";
    end if;
    names := [op(names), objs[i][name]];
  end do;

  num_dof := ComputeDOF(objs, dim);

  if (num_dof < 0) then
    if (nops(hyper_vars) <> -num_dof) then
      candidate_hyp_vars := [];
      for i from 1 to nops(objs) do
        if (objs[i][type] = 'SUPPORT') or (objs[i][type] = 'JOINT') then
          candidate_hyp_vars := [
            op(candidate_hyp_vars),
            op(objs[i][variables])
            ];
        end if;
      end do;
    WARNING(
      "the structure is hyperstatic with %1 overconstrained directions, "
      "please check the structure supports and joints. Also consider defining "
      "the hyperstatic variables by adding 'hyper_vars' property in the "
      "'MakeStructure' method or simply defining '[hyperstatic_variables]' "
      "field in an already existing 'STRUCTURE' object by choosing from the "
      "folloving hyperstatic candidate variables: %2",
      abs(num_dof), candidate_hyp_vars);
    else
      printf("Message (in MakeStructure) "
        "hyperstatic structure detected with %d overconstrained directions\n",
        abs(num_dof) );
    end if;
  elif (num_dof > 0 )then
    error "not enough constraints in the structure";
  else
    printf("Message (in MakeStructure) "
      "isostatic structure detected");
  end if;

  return table({
    type                      = 'STRUCTURE',
    objects                   = objs,
    external_actions          = ext,
    DOF                       = num_dof,
    hyperstatic_variables     = hyper_vars,
    hyperstatic_displacements = hyper_disp,
    dimensions                = dim,
    solved                    = false
    });
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeDOF := proc(
  objs_in::{ # Structure objects
    list({'BEAM', 'ROD', 'SUPPORT', 'JOINT'}),
    set( {'BEAM', 'ROD', 'SUPPORT', 'JOINT'})
  },
  dim::{string}, # Structure dimension ("2D" or "3D")
  $)

  description "Compute the degree of freedom of the input structure objects";

  local dof, i, j, vertex, G, k, objs;
  global EARTH;

  dof  := 0;
  objs := [op(objs_in), EARTH];

  # Built connections graph
  vertex := [];
  printf("Message (in ComputeDOF) checking structure connections... ");
  for i from 1 to nops(objs) do
    vertex := [op(vertex), objs[i][name]];
  end do;
  G := GraphTheory[Graph](vertex);
  for i from 1 to nops(objs) do
    if (objs[i][type] = 'SUPPORT') or (objs[i][type] = 'JOINT') then
      for j from 1 to nops(objs) do
        if member(objs[j][name], objs[i][targets]) then
          GraphTheory[AddEdge](G, {objs[i][name], objs[j][name]});
        end if;
      end do;
    end if;
  end do;

  # Check graph connections
  if GraphTheory[IsConnected](G) then
    printf("\tDONE\n");
  else
    error "unconnected elements detected in the structure";
  end if;

  printf("Message (in ComputeDOF) computing degrees of freedom...");
  for i from 1 to nops(objs) do
    if (objs[i][type] = 'BEAM') then
      if (dim = "2D") then
        dof := dof + 3;
      else
        dof := dof + 6;
      end if;
    elif (objs[i][type] = 'ROD') then
      if (dim = "2D") then
        dof := dof + 3;
      else
        dof := dof + 5;
      end if;
    elif (objs[i][type] = 'JOINT') then
      if (dim = "2D") then
        dof := dof -
          add(objs[i][dof_constr][k], k = [1,2,6]) * (nops(objs[i][targets]) - 1);
      else
        dof := dof -
          add(objs[i][dof_constr][k], k = 1..6) * (nops(objs[i][targets]) - 1);
      end if;
    elif (objs[i][type] = 'SUPPORT') then
      if (dim = "2D") then
        dof := dof -
          add(objs[i][dof_constr][k], k = [1,2,6]) * (nops(objs[i][targets]) - 1);
      else
        dof := dof -
          add(objs[i][dof_constr][k], k = 1..6) * (nops(objs[i][targets]) - 1);
      end if;
    end if;
  end do;
  printf("\tDONE\n");
  printf("Message (in ComputeDOF) display degrees of freedom...", dof);
  printf("\tDOF = %d\n", dof);

  # Display graph
  print( GraphTheory[DrawGraph](G) );

  return dof;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

NE_equilibrium := proc(
  FMQ::{
    list({FORCE,MOMENT,QLOAD}),
    set({FORCE,MOMENT,QLOAD})
  },
  obj::{BEAM,ROD,SUPPORT,JOINT},
  pole,
  dim::string,
  {UP_TO:=obj[length]},
  $)
  
    description "compute the Newton-Euler static equilibrium equations given a set of Forces-Moments-QLoads and the axial coordinate of the pole";
    local eq_T, eq_R, NE_eq, i, comp;
    
    # 2D case
    if dim = "2D" then 
        eq_T := [0,0];
        for i from 1 to nops(FMQ) do
            if FMQ[i][target] = obj[name] then 
                if (FMQ[i][type]=FORCE) then 
                    eq_T := eq_T + FMQ[i][components][1..2] ;
                elif (FMQ[i][type]=QLOAD) then 
                    eq_T := eq_T + map(integrate, FMQ[i][components][1..2](x),x=0..UP_TO);
                end if;
            else
                WARNING("%1 is not applied to %2", FMQ[i], obj);
            end if;
        end do;

        eq_R := [0];
        for i from 1 to nops(FMQ) do
            if FMQ[i][target] = obj[name]  then 
                if (FMQ[i][type]=MOMENT) then 
                    eq_R := eq_R + FMQ[i][components][3] ;
                elif (FMQ[i][type]=FORCE) then 
                    eq_R := eq_R + [FMQ[i][components][2]] *~ (FMQ[i][coordinate] - pole);
                elif (FMQ[i][type]=QLOAD) then 
                    eq_R := eq_R + map(integrate, [FMQ[i][components][2](x)*~(x-pole)], x=0..UP_TO);
                end if;
            else
                WARNING("%1 is not applied to %2", FMQ[i], obj);
            end if;
        end do;
    
    # 3D case
    else
        eq_T := [0,0,0];
        for i from 1 to nops(FMQ) do
            if FMQ[i][target] = obj[name]  then 
                if (FMQ[i][type]=FORCE) then 
                    eq_T := eq_T + FMQ[i][components] ;
                elif (FMQ[i][type]=QLOAD) then 
                    eq_T := eq_T + map(integrate, FMQ[i][components](x),x=0..UP_TO);
                end if;
            else
                WARNING("%1 is not applied to %2", FMQ[i], obj);
            end if;
        end do;

        eq_R := [0,0,0];
        for i from 1 to nops(FMQ) do
            if FMQ[i][target] = obj[name]  then 
                if (FMQ[i][type]=MOMENT) then 
                    eq_R := eq_R + FMQ[i][components] ;
                elif (FMQ[i][type]=FORCE) then 
                    eq_R := eq_R + [0,-FMQ[i][components][3],FMQ[i][components][2]] *~ (FMQ[i][coordinate] - pole);
                elif (FMQ[i][type]=QLOAD) then 
                    eq_R := eq_R + map(integrate, [0,(-FMQ[i][components][3](x)*~(x-pole)),FMQ[i][components][2](x)*~(x-pole)], x=0..UP_TO);
                end if;
            else
                WARNING("%1 is not applied to %2", FMQ[i], obj);
            end if;
        end do;
    end if;   
     
    NE_eq := [op(eq_T),op(eq_R)];
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SolveStructure := proc(struct::STRUCTURE,{compute_intact::boolean := false, compute_disp::boolean := false, compute_stress::boolean := false, shear_contribution::boolean := false, verbose::boolean := false},$)
    local i,g_load,S_obj,S_ext,S_support,S_joint,S_con_forces,vars,eq,sol,obj;
    global _gravity, ground, EARTH;
    description "Solve the Structure
                 INPUTS : Structure object,
                 OUTPUTS: reactions_solutions(list) (supports and joints forces), 
                          reactions FORCES set of FORCE/MOMENT objects,
                          if compute_intact ->update of the internal actions of the structure objects,
                          if compute_disp -> update of the displacement of the structure objects,
                          if compute_stress -> update of the max_VM_stress of the structure objects";
    
    # parsing inputs
    S_obj := {};
    S_ext := {};
    S_support := {};
    S_joint := {};
    S_con_forces := {};
    vars := [];
    for i from 1 to nops(struct[objects]) do
        obj:=struct[objects][i];
        if obj[type] = BEAM or obj[type] = ROD then
            S_obj := {op(S_obj),obj};
        end if;
        if obj[type] = SUPPORT then
            S_support := {op(S_support),obj};
            S_con_forces := {op(S_con_forces),op(obj[forces]),op(obj[moments])};
            vars := [op(vars),op(obj[variables])];
        end if;
        if obj[type] = JOINT then
            S_joint := {op(S_joint),obj};
            S_con_forces := {op(S_con_forces),op(obj[forces]),op(obj[moments])};
            vars := [op(vars),op(obj[variables])];
        end if;
        unassign('obj');
    end do;

    S_ext := struct[external_actions];

    # add gravity distributed load
    if _gravity <> [0,0,0] then
        for i from 1 to nops(S_obj) do
            if S_obj[i][type] = BEAM then
                g_load||(S_obj[i][name]) := MakeQLoad(_gravity*~S_obj[i][cs_area]*~S_obj[i][material][density],S_obj[i],RF=ground);
                S_ext := {op(S_ext),g_load||(S_obj[i][name])};
            end if;
        end do;
    end if; 

    if struct[DOF] = 0 then
        # solve isostatic structure
        printf("Solving the ISOSTATIC structure\n");
        sol:=IsostaticSolver({op(S_obj),op(S_joint),op(S_support)},{op(S_ext),op(S_con_forces)},vars,struct[dimensions],verb=verbose);
        if verbose then 
            printf("solutions:\n");
            print(<sol>);
        end if;
        # update support reactions properties
        printf("updating support reactions fields\n");
        for i from 1 to nops(S_support) do
            S_support[i][support_reactions] := [seq(lhs(S_support[i][support_reactions][j]) = subs(sol,rhs(S_support[i][support_reactions][j])),j=1..nops(S_support[i][support_reactions]))];
        end do;
    elif struct[DOF] < 0 then 
        # solve hyperstatic structure
        if nops(struct[hyperstatic_variables])<>-struct[DOF] then error "mismatch in DOF of the structure, check the hyperstatic variables of the structure and update the structure object"; 
        end if;
        sol:=HyperstaticSolver({op(S_obj),op(S_joint),op(S_support)},{op(S_ext),op(S_con_forces)},vars,struct[hyperstatic_variables],struct[hyperstatic_displacements],struct[dimensions],shear_contrib=shear_contribution,verbosity=verbose);
        if verbose then 
            printf("solutions:\n");
            print(<sol>);
        end if;
        # update support reactions properties
        printf("updating support reactions fields\n");
        for i from 1 to nops(S_support) do
            S_support[i][support_reactions] := [seq(lhs(S_support[i][support_reactions][j]) = subs(sol,rhs(S_support[i][support_reactions][j])),j=1..nops(S_support[i][support_reactions]))];
        end do; 
    end if;

    if compute_intact then
        # computing internal actions
        printf("computing the internal actions\n");
        ComputeInternalActions(S_obj,{op(S_ext),op(S_con_forces)},sol,struct[dimensions],verb = verbose);
    end if;

    # set solved flag to true 
    struct[solved] := true;
    
    return struct;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

HyperstaticSolver := proc(objs::{list({BEAM,ROD,SUPPORT,JOINT}),set({BEAM,ROD})},FMQ::{list({FORCE,MOMENT,QLOAD}),set({FORCE,MOMENT,QLOAD})},vars::list,hyper_vars::list,hyperstatic_displacements::list,dim::string,{shear_contrib::boolean := false,verbosity::boolean := false},$)
    local eq,i, iso_vars, iso_sol, hyper_sol, sol, P, S_obj;
    description "solve hyperstatic structures";
    
    # parse input
    S_obj := {};
    for i from 1 to nops(objs) do
        if objs[i][type] = BEAM or objs[i][type] = ROD then
            S_obj := {op(S_obj),objs[i]};
        end if;
    end do;


    printf("Solving the Hyperstatic variables\n");
    # create a solution in function of hyperstatic variables
    iso_vars := [seq(`if`(member(vars[i],hyper_vars),NULL,vars[i]),i=1..nops(vars))];
    iso_sol:=IsostaticSolver(objs,FMQ,iso_vars,dim,verb=verbosity);

    # compute internal actions
    ComputeInternalActions(S_obj,FMQ,iso_sol,dim,verb = verbosity);
    
    # compute structure internal energy
    P := ComputePotentialEnergy(S_obj,shear_contribution=shear_contrib);

    # hyperstatic equations
    eq:=[];
    for i from 1 to nops(hyper_vars) do
        eq := [op(eq), diff(P,hyper_vars[i]) =  hyperstatic_displacements[i]]
        #eq := [op(eq), coeff(collect(P,hyper_vars[i]),hyper_vars[i]) =  hyperstatic_displacements[i]]
    end do;  
    hyper_sol := op(solve(eq,hyper_vars));

    printf("---DONE---\n");
    sol := [op(hyper_sol),op(subs(hyper_sol,iso_sol))];
    return sol;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputePotentialEnergy := proc(objs::{list({BEAM,ROD}),set({BEAM,ROD})},{shear_contribution:=false},$)
    local i,P;
    description "compute the internal potential energy of the structure";
    
    P := 0;

    for i from 1 to nops(objs) do
        if  member(N,map(lhs,objs[i][internal_actions])) and subs(objs[i][internal_actions](x),N(x))<>0 then
            P := P + integrate(subs(objs[i][internal_actions](x),N(x)^2/(2*objs[i][material][elastic_modulus]*objs[i][cs_area])),x=0..objs[i][length]);
        end if;
        if shear_contribution then
            if member(Ty,map(lhs,objs[i][internal_actions])) and subs(objs[i][internal_actions](x),Ty(x))<>0 then
                P := P + integrate(subs(objs[i][internal_actions](x),objs[i][shear_stiff_factor][1]*Ty(x)^2/(2*objs[i][material][shear_modulus]*objs[i][cs_area])),x=0..objs[i][length]);
            end if;
            if member(Tz,map(lhs,objs[i][internal_actions])) and subs(objs[i][internal_actions](x),Tz(x))<>0 then
                P := P + integrate(subs(objs[i][internal_actions](x),objs[i][shear_stiff_factor][2]*Tz(x)^2/(2*objs[i][material][shear_modulus]*objs[i][cs_area])),x=0..objs[i][length]);
            end if;
        end if;
        if member(Mx,map(lhs,objs[i][internal_actions])) and subs(objs[i][internal_actions](x),Mx(x))<>0 then
            P := P + integrate(subs(objs[i][internal_actions](x),Mx(x)^2/(2*objs[i][material][shear_modulus]*objs[i][cs_inertias][1])),x=0..objs[i][length]);
        end if;
        if member(My,map(lhs,objs[i][internal_actions])) and subs(objs[i][internal_actions](x),My(x))<>0 then
            P := P + integrate(subs(objs[i][internal_actions](x),My(x)^2/(2*objs[i][material][elastic_modulus]*objs[i][cs_inertias][3])),x=0..objs[i][length]);
        end if;
        if member(Mz,map(lhs,objs[i][internal_actions])) and subs(objs[i][internal_actions](x),Mz(x))<>0 then
            P := P + integrate(subs(objs[i][internal_actions](x),Mz(x)^2/(2*objs[i][material][elastic_modulus]*objs[i][cs_inertias][2])),x=0..objs[i][length]);
        end if;
    end do;
    
    return P;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsostaticSolver := proc(objs::{list({BEAM,ROD,SUPPORT,JOINT}),set({BEAM,ROD})},FMQ::{list({FORCE,MOMENT,QLOAD}),set({FORCE,MOMENT,QLOAD})},vars::list,dim::string,{verb::boolean := false},$)
    local eq,i,j,active_FMQ,sol,A,B,rank_eq,vars2;
    description "solve isostatic structures";
    
    # compute equations
    printf("computing the equilibrium equation for the isostatic structure\n");
    eq:=[];
    for i from 1 to nops(objs) do
        active_FMQ := {};
        for j from 1 to nops(FMQ) do
            if FMQ[j][target] = objs[i][name] then
                active_FMQ := {op(active_FMQ), FMQ[j]};
            end if;
        end do;
        eq := [op(eq), op(NE_equilibrium(active_FMQ,objs[i],0,dim))];
        # add joints and supports constraint equations
        if objs[i][type] = SUPPORT or objs[i][type] = JOINT then
            eq := [op(eq),op(objs[i][constraint_loads])];
        end if;
    end do;
    # remove NULL equations
    eq := remove(x -> x = 0 , simplify(eq));
    #remove non used variables
    vars2 := vars;
    for i from 1 to nops(vars) do
        if (has(eq,vars[i])) = false then
            vars2 := remove(x -> x = vars[i], vars2);
            WARNING("%1 was removed from variables because it is not used in the equations",vars[i]);
        end if; 
    end do;
    #equations check;
    A,B := LinearAlgebra[GenerateMatrix](eq,vars2);
    rank_eq := LinearAlgebra[Rank](A);
    if verb then
        printf("structure equilibrium equations:\n");
        print(<op(eq)>);
        printf("structure unknown variables:\n");
        print(vars2);
    end if;
    if rank_eq <> nops(vars2) then error "inconsistent system of equation, got %1 independent equations and %2 variables, check structure supports and joints", rank_eq, nops(vars2); end if;
    # Solve
    printf("solving the structure reaction forces\n");
    sol := simplify(op(solve(eq,vars2)));
    printf("---DONE---\n");
    
    return sol;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeInternalActions := proc(objs::{list({BEAM,ROD}),set({BEAM,ROD})}, FMQ::{list({FORCE,MOMENT,QLOAD}),set({FORCE,MOMENT,QLOAD})},sol::{list,set},dim::string,{verb::boolean := false},$)
    local i, j, active_FMQ, FMQ_sub;
    description "programmatically computation of internal actions for structure objects";

    # substitute structure solution into loads
    FMQ_sub := map2(subs,sol,map(op,FMQ));

    for i from 1 to nops(objs) do 
        # extract active loads
        active_FMQ := {};
        for j from 1 to nops(FMQ_sub) do
            if FMQ_sub[j][target] = objs[i][name] then
                active_FMQ := {op(active_FMQ), FMQ_sub[j]};
            end if;
        end do;
        # compute internal actions
        InternalActions(objs[i], active_FMQ,dim);
    end do;   
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

InternalActions := proc(obj::{BEAM,ROD}, FMQ::{list({FORCE,MOMENT,QLOAD}),set({FORCE,MOMENT,QLOAD})},dim::string,$)
    local i, IA, IAF, IAM, IA_sol, ia_sol, ia_vars, NE_equations, balance,N_,Ty_ ,Tz_ ,Mx_ ,My_ ,Mz_ ;
    description "compute the internal actions of a beam element
                 returns IA as function of 'x' (axial variable)
                 INPUTS : obj  -> object,
                          FMQ  -> list of forces,moments and distributed load";
    
    # printf("Checking actions balance: ");
    # balance := simplify(NE_equilibrium({op(FMQ)},obj,0,dim));
    # if add(abs(balance[i]),i=1..nops(balance))> 0.01 then
    #    error "truss not balanced, equilibrium: %1", balance;
    # else
    #     printf("OK\n");
    # end if:

    # 3D case
    if dim = "3D" then  

    N_ := 0;
    Ty_ := 0;
    Tz_ := 0;
    Mx_ := 0;
    My_ := 0;
    Mz_ := 0;

    # compute internal actions for concentrated loads as effect overlay
    for i from 1 to nops(FMQ) do
        if FMQ[i][type] = FORCE  then
            N_  := `simplify/piecewise`(N_  - piecewise(x>=FMQ[i][coordinate] and x<=obj[length], FMQ[i][components][1]),x) ;
            Ty_ := `simplify/piecewise`(Ty_ + piecewise(x>=FMQ[i][coordinate] and x<=obj[length], FMQ[i][components][2]),x) ;
            Tz_ := `simplify/piecewise`(Tz_ + piecewise(x>=FMQ[i][coordinate] and x<=obj[length], FMQ[i][components][3]),x) ;
            My_ := `simplify/piecewise`(My_ + integrate(piecewise(x>=FMQ[i][coordinate] and x<=obj[length], FMQ[i][components][3]),x=0..x),x) ;
            Mz_ := `simplify/piecewise`(Mz_ + integrate(piecewise(x>=FMQ[i][coordinate] and x<=obj[length], FMQ[i][components][2]),x=0..x),x) ;
        elif FMQ[i][type] = MOMENT then
            Mx_ := `simplify/piecewise`(Mx_ - piecewise(x>=FMQ[i][coordinate] and x<=obj[length], FMQ[i][components][1]),x) ;
            My_ := `simplify/piecewise`(My_ + piecewise(x>=FMQ[i][coordinate] and x<=obj[length], FMQ[i][components][2]),x) ;
            Mz_ := `simplify/piecewise`(Mz_ - piecewise(x>=FMQ[i][coordinate] and x<=obj[length], FMQ[i][components][3]),x) ;
        elif FMQ[i][type] = QLOAD then
            N_  := `simplify/piecewise`(N_  - integrate(FMQ[i][components][1](x),x=0..x),x);
            Ty_ := `simplify/piecewise`(Ty_ + integrate(FMQ[i][components][2](x),x=0..x),x);
            Tz_ := `simplify/piecewise`(Tz_ + integrate(FMQ[i][components][3](x),x=0..x),x);
            My_ := `simplify/piecewise`(My_ + integrate(integrate(FMQ[i][components][3](x),x=0..x),x=0..x),x);
            Mz_ := `simplify/piecewise`(Mz_ + integrate(integrate(FMQ[i][components][2](x),x=0..x),x=0..x),x);
        end if;
    end do;

    IA := [N = unapply(N_,x),
           Ty = unapply(Ty_,x),
           Tz = unapply(Tz_,x),
           Mx = unapply(Mx_,x),
           My = unapply(My_,x),
           Mz = unapply(Mz_,x)
           ];
    
    # 2D case
    elif dim = "2D" then  

    N_ := 0;
    Ty_ := 0;
    Mz_ := 0;
    # compute internal actions for concentrated loads as effect overlay
    for i from 1 to nops(FMQ) do
        if FMQ[i][type] = FORCE  then
            N_  := `simplify/piecewise`(N_  - piecewise(x>=FMQ[i][coordinate] and x<=obj[length], FMQ[i][components][1]),x) ;
            Ty_ := `simplify/piecewise`(Ty_ + piecewise(x>=FMQ[i][coordinate] and x<=obj[length], FMQ[i][components][2]),x) ;
            Mz_ := `simplify/piecewise`(Mz_ + integrate(piecewise(x>=FMQ[i][coordinate] and x<=obj[length], FMQ[i][components][2]),x=0..x),x) ;
        elif FMQ[i][type] = MOMENT then
            Mz_ := `simplify/piecewise`(Mz_ - piecewise(x>=FMQ[i][coordinate] and x<=obj[length], FMQ[i][components][3]),x) ;
        elif FMQ[i][type] = QLOAD then       
            N_  := `simplify/piecewise`(N_  - integrate(FMQ[i][components][1](x),x=0..x),x);
            Ty_ := `simplify/piecewise`(Ty_ + integrate(FMQ[i][components][2](x),x=0..x),x);
            Mz_ := `simplify/piecewise`(Mz_ + integrate(integrate(FMQ[i][components][2](x),x=0..x),x=0..x),x);        end if;
    end do;

    IA := [N = unapply(N_,x),
           Ty = unapply(Ty_,x),
           Mz = unapply(Mz_,x)
           ];

    end if;

    if obj[type]=ROD then
        IA := [IA[1]];
    end if;

    printf("updating %s %s's internal actions\n",obj[type], obj[name]);
    obj[internal_actions] := IA;
    printf("---DONE---\n");

    return ;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# That's all folks!
