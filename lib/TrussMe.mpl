# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                     _____                   __  __                          #
#                    |_   _| __ _   _ ___ ___|  \/  | ___                     #
#                      | || '__| | | / __/ __| |\/| |/ _ \                    #
#                      | || |  | |_| \__ \__ \ |  | |  __/                    #
#                      |_||_|   \__,_|___/___/_|  |_|\___|                    #
#                A Maple Library for Truss Elements Structures                #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Authors:
#  Matteo Larcher (University of Trento)
#  Davide Stocco (University of Trento)
#
# License: BSD 3-Clause License
#
# This is a module for the 'TrussMe' (A Maple Library for Truss Elements
# Structures).

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

TrussMe := module()

# TODO: Hide module content with dummy procedures as procname := proc(inputs); _procname(inputs); end proc;
#       So the source code is not visible by showstat and showsource commands.
export  `union`,
        SetModuleOptions,
        PrintStartProc,
        PrintEndProc,
        IsEarth,
        SetGravity,
        GetGravity,
        Show,
        AssignData,
        UnAssignData,
        Rotate,
        Translate,
        Project,
        InverseFrame,
        IsFrame,
        IsPoint,
        IsVector,
        Origin,
        Uvec,
        UvecX,
        UvecY,
        UvecZ,
        Norm2,
        MakeMaterial,
        IsMaterial,
        MakeBeam,
        MakeBeamPoints,
        IsBeam,
        MakeRod,
        MakeRodPoints,
        IsRod,
        MakeRigidBody,
        IsRigidBody,
        MakeJoint,
        IsJoint,
        MakeSupport,
        IsSupport,
        IsCompliantSupport,
        IsCompliantJoint,
        MakeForce,
        IsForce,
        MakeMoment,
        IsMoment,
        MakeQForce,
        IsQForce,
        MakeQMoment,
        IsQMoment,
        MakeStructure,
        IsStructure,
        SolveStructure,
        PlotStructure,
        PlotDeformedStructure,
        CleanJoint,
        CleanSupport,
        CleanRod,
        CleanBeam,
        CleanStructure,
        DrawStructureGraph,
        DrawStructureSparseMatrix,
        ComputePunctualDisplacement;

global  ground;

local   ModuleLoad,
        ModuleUnload,
        earth,
        gravity,
        GetNames,
        GetObjByName,
        GetObjsByType,
        CopyStructure,
        Simplify,
        Subs,
        Diff,
        ComputeDOF,
        NewtonEuler,
        HyperstaticSolver,
        IsostaticSolver,
        LinearSolver,
        ComputeInternalActions,
        ComputePotentialEnergy,
        ComputeDisplacements,
        InternalActions,
        InitTrussMe,
        TypeRegister,
        Protect,
        ComputeSpringDisplacement,
        ComputeSpringEnergy,
        ComputeSupportDisplacements,
        ComputeJointDisplacements,
        ComputeObjectFrameDisplacements,
        ObjectColor,
        PlotRigidBody,
        PlotDeformedRigidBody,
        PlotBeam,
        PlotDeformedBeam,
        PlotRod,
        PlotDeformedRod,
        PlotJoint,
        PlotDeformedJoint,
        PlotSupport,
        PlotDeformedSupport,
        IsInsideJoint,
        IsInsideSupport,
        IsInsideBeam,
        IsInsideRod,
        IsInsideStructure,
        lib_base_path,
        verbose_mode,
        suppress_warnings,
        time_limit_simplify,
        print_indent,
        print_increment,
        ListPadding,
        keep_veiled,
        Beam_color,
        Rod_color,
        RigidBody_color,
        CompliantSupport_color,
        Support_color,
        CompliantJoint_color,
        Joint_color,
        Earth_color,
        StoredData,
        veiling_label;

option  package,
        load   = ModuleLoad,
        unload = ModuleUnload;

description "A Maple Library for Truss Elements Structures.";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   __  __           _       _      _                    _
#  |  \/  | ___   __| |_   _| | ___| |    ___   __ _  __| |
#  | |\/| |/ _ \ / _` | | | | |/ _ \ |   / _ \ / _` |/ _` |
#  | |  | | (_) | (_| | |_| | |  __/ |__| (_) | (_| | (_| |
#  |_|  |_|\___/ \__,_|\__,_|_|\___|_____\___/ \__,_|\__,_|

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ModuleLoad := proc()

  description "Module 'TrussMe' load procedure.";

  local i;

  # Display module init message
  printf(
    "'TrussMe' module version beta-0.0, BSD 3-Clause License - "
    "Copyright (c) 2023, Matteo Larcher and Davide Stocco.\n"
  );

  # Library path
  lib_base_path := null;
  for i in [libname] do
    if (StringTools[Search]("TrussMe", i) <> 0) then
      lib_base_path := i;
    end if;
  end do;
  if (lib_base_path = null) then
    error "Cannot find 'TrussMe' library.";
  end if;

  # Register types
  TypeRegister();

  # Initialize the module variables
  InitTrussMe();

  # Protect Module Keywords
  Protect();

  return NULL;
end proc: # ModuleLoad

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   __  __           _       _      _   _       _                 _
#  |  \/  | ___   __| |_   _| | ___| | | |_ __ | | ___   __ _  __| |
#  | |\/| |/ _ \ / _` | | | | |/ _ \ | | | '_ \| |/ _ \ / _` |/ _` |
#  | |  | | (_) | (_| | |_| | |  __/ |_| | | | | | (_) | (_| | (_| |
#  |_|  |_|\___/ \__,_|\__,_|_|\___|\___/|_| |_|_|\___/ \__,_|\__,_|

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ModuleUnload := proc()

  description "Module 'TrussMe' unload procedure.";

  #printf("Unloading 'TrussMe'\n");
  verbose_mode      := 1;
  suppress_warnings := false;
  UnAssignData();

end proc: # ModuleUnload

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   _____                 ____            _     _
#  |_   _|   _ _ __   ___|  _ \ ___  __ _(_)___| |_ ___ _ __
#    | || | | | '_ \ / _ \ |_) / _ \/ _` | / __| __/ _ \ '__|
#    | || |_| | |_) |  __/  _ <  __/ (_| | \__ \ ||  __/ |
#    |_| \__, | .__/ \___|_| \_\___|\__, |_|___/\__\___|_|
#        |___/|_|                   |___/

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

TypeRegister := proc()

  description "Register 'TrussMe' module types.";

  # Register types
  TypeTools[AddType](EARTH, IsEarth);
  TypeTools[AddType](FRAME, IsFrame);
  TypeTools[AddType](POINT, IsPoint);
  TypeTools[AddType](VECTOR, IsVector);
  TypeTools[AddType](BEAM, IsBeam);
  TypeTools[AddType](ROD, IsRod);
  TypeTools[AddType](RIGID_BODY, IsRigidBody);
  TypeTools[AddType](FORCE, IsForce);
  TypeTools[AddType](MOMENT, IsMoment);
  TypeTools[AddType](QFORCE, IsQForce);
  TypeTools[AddType](QMOMENT, IsQMoment);
  TypeTools[AddType](SUPPORT, IsSupport);
  TypeTools[AddType](JOINT, IsJoint);
  TypeTools[AddType](MATERIAL, IsMaterial);
  TypeTools[AddType](STRUCTURE, IsStructure);

end proc: # TypeRegister

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   ___       _ _  _____                   __  __
#  |_ _|_ __ (_) ||_   _| __ _   _ ___ ___|  \/  | ___
#   | || '_ \| | __|| || '__| | | / __/ __| |\/| |/ _ \
#   | || | | | | |_ | || |  | |_| \__ \__ \ |  | |  __/
#  |___|_| |_|_|\__||_||_|   \__,_|___/___/_|  |_|\___|

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

InitTrussMe := proc()

  description "Initialize 'TrussMe' module internal variables.";

  # Define Module Variables
  ground := <<1, 0, 0, 0>|
             <0, 1, 0, 0>|
             <0, 0, 1, 0>|
             <0, 0, 0, 1>>;

  gravity := [0, 0, 0];

  verbose_mode           := 1;
  suppress_warnings      := false;
  time_limit_simplify    := 5;
  print_indent           := 0;
  print_increment        := 4;
  Beam_color             := "SteelBlue";
  Rod_color              := "Niagara DarkOrchid";
  RigidBody_color        := "Indigo";
  CompliantSupport_color := "DarkGreen";
  Support_color          := "DarkOrange";
  CompliantJoint_color   := "LightSalmon";
  Joint_color            := "MediumSeaGreen";
  Earth_color            := "Firebrick";
  StoredData             := [];
  veiling_label          := '_V';

  earth := table({
    "type"             = EARTH,
    "name"             = "earth",
    "length"           = 0,
    "frame"            = ground,
    "admissible_loads" = [1, 1, 1, 1, 1, 1]
    }):

end proc: # InitTrussMe

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   ____            _            _
#  |  _ \ _ __ ___ | |_ ___  ___| |_
#  | |_) | '__/ _ \| __/ _ \/ __| __|
#  |  __/| | | (_) | ||  __/ (__| |_
#  |_|   |_|  \___/ \__\___|\___|\__|

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Protect := proc()

  description "Protect 'TrussMe' module internal variables.";

  # Protect Module Global Variables
  protect(
    'ground'
  );

  # Protect the types
  protect(
    'FRAME',
    'EARTH',
    'BEAM',
    'ROD',
    'FORCE',
    'MOMENT',
    'QFORCE',
    'QMOMENT',
    'SUPPORT',
    'JOINT',
    'MATERIAL',
    'STRUCTURE'
  );

end proc: # Protect

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#    ___                      _     _
#   / _ \__   _____ _ __ _ __(_) __| | ___  ___
#  | | | \ \ / / _ \ '__| '__| |/ _` |/ _ \/ __|
#  | |_| |\ V /  __/ |  | |  | | (_| |  __/\__ \
#   \___/  \_/ \___|_|  |_|  |_|\__,_|\___||___/
#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

`union` := proc(
  A::{list, set}, # Object A to be united
  B::{list, set}, # Object B to be united
  $)::{list, set};

  option overload;

  description "Extension of union operator to list objects <A> and <B>.";

  local out;
  PrintStartProc(procname);

  if type(A, 'set') and type(B, 'set') then
    out := {op(A), op(B)};
  else
    out := [op(A), op(B)];
  end if;

  PrintEndProc(procname);
  return out;
end proc: # union

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   _____                 _   _
#  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___
#  | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
#  |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
#  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SetModuleOptions := proc(
  {
    verbosity::{integer, nothing}        := NULL, # Verbose mode
    disable_warnings::{boolean, nothing} := NULL, # Suppress warnings
    time_limit::{constant, nothing}      := NULL  # Time limit for simplify operations
  },
  $)::{nothing};

  description "Set the module options: <verbose_mode>::integer = [0, 1, 2], "
    "<disable_warnings>::boolean = [true, false] and <time_limit>::constant"
    " = [0, +inf].";

  # Verbosity
  if (verbosity <> NULL) then
    if (verbosity < 0) or (verbosity > 2) then
      error "invalid verbose mode detected";
    else
      verbose_mode := verbosity;
    end if;
  end if;

  # Suppress warnings
  if (disable_warnings <> NULL) then
    suppress_warnings := disable_warnings;
  end if;

  # Time limit
  if (time_limit <> NULL) then
    if (time_limit < 0) then
      error "invalid time limit detected";
    else
      time_limit_simplify := time_limit;
    end if;
  end if;

  return NULL;
end proc: # SetModuleOptions

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PrintStartProc := proc(
  proc_name::{procedure}, # Procedure name
  $)::{nothing};

  description "Print the start message of a procedure with name <proc_name>.";

  # Increase printf indentation
  print_indent := print_indent + print_increment;

  # Show start message
  if (verbose_mode > 1) then
    printf("%*sStart '%s' procedure...\n", print_indent, "|   ", proc_name);
  end if;
end proc: # PrintStartProc


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PrintEndProc := proc(
  proc_name::{procedure}, # Procedure name
  $)::{nothing};

  description "Print the end message of a procedure with name <proc_name>.";

  # Show end message
  if (verbose_mode > 1) then
    printf("%*sEnd   '%s' procedure\n", print_indent, "|   ", proc_name);
  end if;

  # Increase printf indentation
  print_indent := print_indent - print_increment;
end proc: # PrintEndProc

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsEarth := proc(
  obj::{anything}, # Object to be tested
  $)::{boolean};

  description "Test if an object <obj> is the EARTH object.";

  PrintStartProc(procname);
  PrintEndProc(procname);
  return evalb(obj["type"] = EARTH);
end proc: # IsEarth

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SetGravity := proc(
  obj::{VECTOR}, # Gravity vector
  $)::{nothing};

  description "Set gravity vector with [x, y, z] components of <obj>.";

  PrintStartProc(procname);
  # Set gravity local variable
  if (nops(obj) = 3) then
    gravity := obj;
  else
    error "invalid gravity vector detected";
  end if;
  PrintEndProc(procname);
  return NULL;
end proc: # SetGravity

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

GetGravity := proc( $ )::{VECTOR};

  description "Get gravity vector.";

  PrintStartProc(procname);
  PrintEndProc(procname);
  return gravity;
end proc: # GetGravity

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Norm2 := proc(
  obj::{list, vector}, # Vector for which the norm is computed
  $)::{algebraic};

  description "Compute the norm of a vector <obj>.";

  local out;
  PrintStartProc(procname);

  out := sqrt(add(x, x in obj^~2));

  PrintEndProc(procname);
  return out;
end proc: # Norm2

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ListPadding := proc(
  lst::{list, Vector, algebraic}, # List to be padded
  n::{integer},                   # Number of elements of the final list
  value::{algebraic} := 0,        # Value to be used for padding
  $)::{list, Vector};

  description "Pad a list <list> with <value> to have <n> elements.";

  local i, out;
  PrintStartProc(procname);

  if type(lst, algebraic) then
    out := [lst];
  elif type(lst, Vector) then
    out := [op(convert(lst, list))];
  else
    out := lst;
  end if;

  if (nops(out) < n) then
    out := out union [seq(value, i = (1..n-nops(out)))];
  elif (nops(out) > n) then
    out := out[1..n];
  end if;

  if type(lst, Vector) then
    out := convert(out, Vector);
  end if;

  PrintEndProc(procname);
  return out;
end proc: # ListPadding

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Show := proc(
  tab::{table}, # Table to be shown
  $)::{nothing};

  description "Show the content of a table <tab>.";

  PrintStartProc(procname);
  print(tab = tab["type"](op(op(tab))));
  PrintEndProc(procname);
end proc: # Show

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

GetNames := proc(
  objs::{ # Structural elements
    list({MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH}),
    set( {MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH})
  }, $)::{list({string}), set({string})};

  description "Get names of a list/set of objects <objs>.";

  local out;
  PrintStartProc(procname);

  if type(objs, 'set') then
    out := {seq(objs[i]["name"], i = 1..nops(objs))};
  else
    out := [seq(objs[i]["name"], i = 1..nops(objs))];
  end if;

  PrintEndProc(procname);
  return out;
end proc: # GetNames

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

GetObjByName := proc(
  name::{string}, # Name of the object
  objs::{         # Structural elements
    list({MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH}),
    set( {MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH})
  }, $)::{MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH};

  description "Get object which name field is <name> from a list/set of objects "
    "<objs>.";

  local out, obj;
  PrintStartProc(procname);

  out := NULL;

  for obj in objs do
    if (obj["name"] = name) then
      out := eval(obj); # Do not remove eval
      break;
    end if;
  end do;

  PrintEndProc(procname);
  return out;
end proc: # GetObjByName

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

GetObjsByType := proc(
  types::{list(symbol), set(symbol)}, # List of types to be selected from <objs>
  objs::{list, set},                  # Structural elements
  $)::{list};

  description "Get objects which type field is in <types> from a list/set of "
    "objects <objs>.";

  local out, obj;
  PrintStartProc(procname);

  out := [];
  for obj in objs do
    if (obj::convert(types, set)) then
      out := out union [eval(obj)]; # Do not remove eval
    end if;
  end do;

  PrintEndProc(procname);
  return out;
end proc: # GetObjsByType

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Simplify := proc(
  obj::{anything},         # Expression to be simplified
  opt::{anything} := NULL, # Options
  $)::{anything};

  description "Simplify an algebraic expression <obj>.";

  local out, time_limit;
  PrintStartProc(procname);

  time_limit := `if`(procname::indexed, op(procname), time_limit_simplify);

  try
    timelimit(time_limit, simplify(obj, opt));
    out := %;
  catch :
    WARNING("Time limit of %1s exceeded for simplify operation, raw solutions "
      "is returned. The input <time_limit> can be modified by setting it in the "
      "SetModuleOptions procedure.", time_limit
    );
    out := obj;
  end try:

  PrintEndProc(procname);
  return eval(out);
end proc: # Simplify

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  AssignData := proc(
    x::{list, set}, # The list to be assigned
    $)::{nothing};

    description "Assign the list <x> to the local variable <StoredData>.";

    StoredData := x;
    return NULL;
  end proc: # AssignData

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  UnAssignData := proc(
    $)::{nothing};

    description "Unassign the local variable <StoredData>.";

    StoredData := [];
    return NULL;
  end proc: # UnAssignData

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Subs := proc()::{anything};

  description "Perform subs command neglecting sub-lists and sub-sets from the "
    "substitution list.";

  local x, y, out;
  PrintStartProc(procname);

  map(x -> map(remove, y-> type(y, {list, set}), x), [_passed[1..-2]]);
  out := subs(op(%), _passed[-1]);

  PrintEndProc(procname);
  return out;
end proc; # Subs

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Diff := proc({veils := NULL})::anything;

  description "Perform diff command on veiled expressions given the veiling list "
    "<veils>.";

  local out, subs_diff, d_vars, v2f, f2v, last;
  PrintStartProc(procname);

  if (veils = NULL) then
    veils := [];
    last := -1
  else
    last := -2;
  end if;

  # Get the variables to be differentiated
  d_vars := _passed[2..last];

  # Veil to functions substitution list
  v2f := map(x -> lhs(x) =~ lhs(x)(d_vars), veils);

  # Function to veils substitution list
  f2v := rhs~(v2f) =~ lhs~(v2f);

  subs(v2f, veils);
  diff(lhs~(%), d_vars) =~ Simplify(diff(rhs~(%), d_vars));
  subs_diff := lhs~(%) =~ Simplify(subs(op(ListTools[Reverse](%)),rhs~(%))):

  # Compute the derivative of the veiled expression
  subs(subs_diff, diff(subs(v2f, _passed[1]), d_vars));

  # Substitute back the veils
  out := subs(f2v, %);

  PrintEndProc(procname);
  return out;
end proc; # Diff

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

InverseFrame := proc(
  RF::{FRAME}, # Reference frame (affine transformation) to be inverted
  $)::{FRAME};

  description "Inverse transformation matrix of an affine transformation <RF>.";

  local out;
  PrintStartProc(procname);

  LinearAlgebra[Transpose](RF[1..3, 1..3]);
  out := <<% | -% . RF[1..3, 4]>,
          <0 | 0 | 0 | 1>>;

  PrintEndProc(procname);
  return out;
end proc: # InverseFrame

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsFrame := proc(
  obj::{anything}, # Object to be checked
  $)::{boolean};

  description "Check if the input object <obj> is a FRAME object.";

  local out;
  PrintStartProc(procname);

  if (type(obj, 'Matrix')) and
     (LinearAlgebra[RowDimension](obj) = 4) and
     (LinearAlgebra[ColumnDimension](obj) = 4) then
    out := true;
  else
    out := false;
  end if;

  PrintEndProc(procname);
  return out;
end proc: # IsFrame

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Rotate := proc(
  axis::{symbol, string}, # Rotation axis
  angle::{algebraic},     # Rotation angle (rad)
  $)::{FRAME};

  description "Transformation matrix corresponding to the rotation <angle> "
    "around the given <axis>";

  local out;
  PrintStartProc(procname);

  if (axis = 'X') or (axis = "X") then
    out := <<1, 0,           0,          0>|
            <0, cos(angle),  sin(angle), 0>|
            <0, -sin(angle), cos(angle), 0>|
            <0, 0,           0,          1>>;
  elif (axis = 'Y') or (axis = "Y") then
    out := <<cos(angle), 0, -sin(angle), 0>|
            <0,          1, 0,           0>|
            <sin(angle), 0, cos(angle),  0>|
            <0,          0, 0,           1>>;
  elif (axis = 'Z') or (axis = "Z") then
    out := <<cos(angle),  sin(angle), 0, 0>|
            <-sin(angle), cos(angle), 0, 0>|
            <0,           0,          1, 0>|
            <0,           0,          0, 1>>;
  else
    error "invalid axis detected";
  end if;

  PrintEndProc(procname);
  return out;
end proc: # Rotate

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Translate := proc(
  x::{algebraic}, # X-axis translation component
  y::{algebraic}, # Y-axis translation component
  z::{algebraic}, # Z-axis translation component
  $)::{FRAME};

  description "Transformation matrix corresponding to the translation <x,y,z>.";

  local out;
  PrintStartProc(procname);

  out := <<1, 0, 0, 0>|
          <0, 1, 0, 0>|
          <0, 0, 1, 0>|
          <x, y, z, 1>>;

  PrintEndProc(procname);
  return out;
end proc: # Translate

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Origin := proc(
  RF::{FRAME}, # Reference frame
  $)::{vector};

  description "Extract the origin of the reference frame <RF>.";

  local out;
  PrintStartProc(procname);

  out := <RF[1,4], RF[2,4], RF[3,4], 1>;

  PrintEndProc(procname);
  return out;
end proc: # Translate

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsPoint := proc(
  obj::{anything}, # Object to be checked
  $)::{boolean};

  description "Check if the input object <obj> is a POINT object.";

  local out;
  PrintStartProc(procname);

  if (type(obj, 'list')) and
     (nops(obj) = 3) then
    out := true;
  else
    out := false;
  end if;

  PrintEndProc(procname);
  return out;
end proc: # IsPoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsVector := proc(
  obj::{anything}, # Object to be checked
  $)::{boolean};

  description "Check if the input object <obj> is a VECTOR object.";

  local out;
  PrintStartProc(procname);

  if (type(obj, 'list')) and
     (nops(obj) = 3) then
    out := true;
  else
    out := false;
  end if;

  PrintEndProc(procname);
  return out;
end proc: # IsVector

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Uvec := proc(
  axis::{symbol},        # Axis of the unit vector
  RF::{FRAME} := ground, # Reference frame
  $)::{Vector};

  description "Extract the unit vector of the reference frame <RF> along the "
    "given <axis>.";

  local out;
  PrintStartProc(procname);

  if (axis = 'X') then
    out := <RF[1,1], RF[2,1], RF[3,1], 0>;
  elif (axis = 'Y') then
    out := <RF[1,2], RF[2,2], RF[3,2], 0>;
  elif (axis = 'Z') then
    out := <RF[1,3], RF[2,3], RF[3,3], 0>;
  else
    error "invalid axis detected";
  end if;

  PrintEndProc(procname);
  return out;
end proc: # Uvec

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

UvecX := proc(
  RF::{FRAME} := ground, # Reference frame
  $)::{Vector};

  description "Extract the x-axis unit vector of the reference frame <RF>.";

  local out;
  PrintStartProc(procname);

  out := <RF[1,1], RF[2,1], RF[3,1], 0>;

  PrintEndProc(procname);
  return out;
end proc: # UvecX

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

UvecY := proc(
  RF::{FRAME} := ground, # Reference frame
  $)::{Vector};

  description "Extract the y-axis unit vector of the reference frame <RF>.";

  local out;
  PrintStartProc(procname);

  out := <RF[1,2], RF[2,2], RF[3,2], 0>;

  PrintEndProc(procname);
  return out;
end proc: # UvecY

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

UvecZ := proc(
  RF::{FRAME} := ground, # Reference frame
  $)::{Vector};

  description "Extract the z-axis unit vector of the reference frame <RF>.";

  local out;
  PrintStartProc(procname);

  out := <RF[1,3], RF[2,3], RF[3,3], 0>;

  PrintEndProc(procname);
  return out;
end proc: # UvecZ

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Project := proc(
  x::{list, Vector}, # Object to be projected
  RF_ini::{FRAME},   # Reference frame from which the object is expressed
  RF_end::{FRAME},   # Reference frame to which the object will be expressed
  $)::{list, Vector};

  description "Project <x,y,z>, or vector <x,y,z,0>, or point <x,y,z,1> from "
    "reference frame <RF_ini> to reference frame <RF_end>.";

  local x_tmp, out;

  PrintStartProc(procname);

  # Pad input vector with 0 if its length is 3
  if not (nops(x) = 3) and
     not (nops(x) = 4) then
    error "invalid input vector/point <x> detected";
  end if;
  x_tmp := ListPadding(convert(x, Vector), 4);

  try
    # try to compare RF_end and RF_ini
    # FIXME: problems with floats (floats not handled error)
    map(evalb, evala(simplify(RF_end)) =~ evala(simplify(RF_ini)));
  catch:
    map(evalb, (RF_end) =~ (RF_ini));
  end try;
  if has(%, false) then
    InverseFrame(RF_end).RF_ini.x_tmp;
    out := Simplify([seq(%[i], i = 1..nops(x))]);
  else
    out := x_tmp[1..nops(x)];
  end if;

  out := convert(out, whattype(x));

  PrintEndProc(procname);
  return out;
end proc: # Project

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeMaterial := proc({
    name::{string}               := "DeafultSteel", # Name of the material
    elastic_modulus::{algebraic} := 210.0E+09,      # Elastic modulus (Pa)
    poisson_ratio::{algebraic}   := 0.3,            # Poisson ratio (-)
    shear_modulus::{algebraic}   := elastic_modulus/(2*(1+poisson_ratio)),
                                                    # Shear modulus (Pa)
    density::{algebraic}         := 7.4E+03         # Density (kg/m^3)
  }, $)::{MATERIAL};

  description "Define a MATERIAL object with inputs: name of the material, "
    "elastic modulus <elastic_modulus> (default = 210.0E9 Pa), Poisson ratio "
    "<poisson_ratio> (default = 0.3), shear modulus <shear_modulus> (default "
    "= E/(2*(1+nu))), density <density> (default = 7.4E3 kg/m^3).";

  local out;
  PrintStartProc(procname);

  out := table({
    "type"            = MATERIAL,
    "name"            = name,
    "elastic_modulus" = elastic_modulus,
    "poisson_ratio"   = poisson_ratio,
    "shear_modulus"   = shear_modulus,
    "density"         = density
    });

  PrintEndProc(procname);
  return op(out);
end proc: # DefineMaterial

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsMaterial := proc(
  obj::{anything}, # Object to be checked
  $)::{boolean};

  description "Check if the input object <obj> is a MATERIAL object.";

  PrintStartProc(procname);
  PrintEndProc(procname);
  return evalb(obj["type"] = MATERIAL);
end proc: # IsMaterial

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeForce := proc(
  components::{VECTOR},                                # Force components in RF
  coords::{algebraic, list(algebraic)},                # Application coordinates in object frame
  obj::{BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH}, # Target object
  RF::{FRAME} := ground,                               # Reference frame
  $)::{FORCE};

  description "Define a FORCE object with inputs: force components <components>, "
    "force application axial coordinate <coords>, target object <obj>, and "
    "optional reference frame <RF> in which the force is defined (default = "
    "ground).";

  local proj_components, admissible_components, out;

  PrintStartProc(procname);

  proj_components       := Project(components, RF, obj["frame"]);
  admissible_components := convert(proj_components .~ <obj["admissible_loads"][1..3]>, list);

  # Check input arguments
  if (proj_components <> admissible_components) and (not suppress_warnings) then
  ["x_comp", "y_comp", "z_comp"] =~
    convert(proj_components .~ <eval(map((x->evalb(x = 0)), obj["admissible_loads"][1..3]), [true = 1, false = 0])>, list);
    WARNING("Force components are not admissible for the target object. The "
      "following components will be ignored: %1", remove(x-> rhs(x) = 0, %));
  end if;

  if IsSupport(obj) or IsJoint(obj) then
    if (ListPadding(coords, 3) <> [0,0,0]) then
      error "only null axial coordinate is accepted for SUPPORT and JOINT "
        "objects";
    end if;
  end if;

  out := table({
    "type"       = FORCE,
    "components" = admissible_components,
    "coordinate" = ListPadding(coords, 3),
    "target"     = obj["name"]
    });

  PrintEndProc(procname);
  return op(out);
end proc: # MakeForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsForce := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a FORCE object";

  PrintStartProc(procname);
  PrintEndProc(procname);
  return evalb(obj["type"] = FORCE);
end proc: # IsForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeMoment := proc(
  components::{VECTOR},                           # Moment components
  coords::{algebraic, list(algebraic)},           # Application coordinates in object frame
  obj::{BEAM, RIGID_BODY, SUPPORT, JOINT, EARTH}, # Target object
  RF::{FRAME} := ground,                          # Reference frame
  $)::{MOMENT};

  description "Define a MOMENT object with inputs: moment components <components>, "
    "moment application axial coordinate <coords>, target object <obj>, and "
    "optional reference frame <RF> in which the moment is  defined (default = "
    "ground).";

  local proj_components, admissible_components, out;
  PrintStartProc(procname);

  proj_components       := Project(components, RF, obj["frame"]);
  admissible_components := convert(proj_components .~ <obj["admissible_loads"][4..6]>, list);

  # Check input arguments
  if (proj_components <> admissible_components) and (not suppress_warnings) then
    ["x_comp", "y_comp", "z_comp"] =~
      convert(proj_components .~ <eval(map((x->evalb(x = 0)), obj["admissible_loads"][4..6]), [true = 1, false = 0])>, list);
    WARNING("Moment components are not admissible for the target object. The "
      "following components will be ignored: %1", remove(x-> rhs(x) = 0, %));
  end if;

  if IsSupport(obj) or IsJoint(obj) then
    if (ListPadding(coords, 3) <> [0,0,0]) then
      error "only null axial coordinate is accepted for SUPPORT and JOINT "
        "objects";
    end if;
  end if;

  out := table({
    "type"       = MOMENT,
    "components" = admissible_components,
    "coordinate" = ListPadding(coords, 3),
    "target"     = obj["name"]
    });

  PrintEndProc(procname);
  return op(out);
end proc: # MakeMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsMoment := proc(
  obj::{anything}, # Object to be checked
  $)::{boolean};

  description "Check if the object <obj> is a MOMENT object.";

  PrintStartProc(procname);
  PrintEndProc(procname);
  return evalb(obj["type"] = MOMENT);
end proc: # IsMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeQForce := proc(
  components::{procedure, list(algebraic)}, # Distributed load components
  obj::{BEAM, ROD},                         # Target object
  RF::{FRAME} := ground,                    # Reference frame
  {
    ell_min::{algebraic} := 0,              # Initial axial coordinate
    ell_max::{algebraic} := obj["length"]   # Final axial coordinate
  }, $)::{QFORCE};

  description "Define a QFORCE object with inputs: distributed load target "
    "object components <components>, target object <obj>, optional reference "
    "frame <RF> in which the load components are defined (default = ground), "
    "and optional initial <ell_min> and final <ell_max> application points "
    "(axial coordinates).";

  local proj_components, x, out;
  PrintStartProc(procname);

  if type(components, procedure) then
    proj_components := unapply(Project(components(x), RF, obj["frame"]), x);
  else
    proj_components := (x) -> piecewise((ell_min <= x) and (x <= ell_max), Project(components, RF, obj["frame"]), 0);
  end if;

  if IsRod(obj) then
    if (proj_components(x)[2] <> 0) or (proj_components(x)[3] <> 0) then
      error "only axial loads are accepted in ROD objects"
    end if;
  end if;

  out := table({
    "type"        = QFORCE,
    "components"  = proj_components,
    "target"      = obj["name"]
    });

  PrintEndProc(procname);
  return op(out);
end proc: # MakeQForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsQForce := proc(
  obj::{anything}, # Object to be checked
  $)::{boolean};

  description "Check if the object <obj> is a QFORCE object.";

  PrintStartProc(procname);
  PrintEndProc(procname);
  return evalb(obj["type"] = QFORCE);
end proc: # IsQForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeQMoment := proc(
  components::{procedure, list(algebraic)}, # Distributed load components
  obj::{BEAM},                              # Target object
  RF::{FRAME} := ground,                    # Reference frame in which the moment is defined
  {
    ell_min::{algebraic} := 0,              # Initial application point (axial coordinate)
    ell_max::{algebraic} := obj["length"]   # Final application point (axial coordinate)
  }, $)::{QMOMENT};

  description "Define a QMOMENT object with inputs: distributed torque target "
    "object components <components>, target object <obj>, optional reference "
    "frame <RF> in which the load components are defined (default = ground), "
    "and optional initial <ell_min> and final <ell_max> application points "
    "(axial coordinates).";

  local proj_components, x, out;
  PrintStartProc(procname);

  if type(components, procedure) then
    proj_components := unapply(Project(components(x), RF, obj["frame"]),x);
  else
    proj_components := (x) -> piecewise((ell_min <= x) and (x <= ell_max), Project(components, RF, obj["frame"]), 0);
  end if;

  out := table({
    "type"        = QMOMENT,
    "components"  = proj_components,
    "target"      = obj["name"]
    });

  PrintEndProc(procname);
  return op(out);
end proc: # MakeQMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsQMoment := proc(
  obj::{anything}, # Object to be checked
  $)::{boolean};

  description "Check if the object <obj> is a QMOMENT object";

  PrintStartProc(procname);
  PrintEndProc(procname);
  return evalb(obj["type"] = QMOMENT);
end proc: # IsQMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeSupport := proc(
  name::{string},                                        # Support name
  constrained_dof::{list},                               # Constrained degree of freedom
  objs::{list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})}, # Target objects
  coords::{list},                                        # Support locations
  RF::{FRAME} := ground,                                 # Reference frame of the support
  {
    stiffness::{procedure,list(algebraic)} := [ # Stiffness components (default = infinite)
      infinity, infinity, infinity,
      infinity, infinity, infinity
    ] *~ constrained_dof
  }, $)::{SUPPORT};

  description "Make a SUPPORT object with inputs: support name <name>, constrained "
    "degrees of freedom <constrained_dof>, target objects <objs>, support locations "
    "<coords>, and optional reference frame <RF> in which the support is defined "
    "(default = ground). The optional input <stiffness> is a list of stiffness "
    "components (default = infinite) in the order: [ktx, kty, ktz, krx, kry, krz].";

  local S, J_tmp, i, j, sr_F_names, sr_F_values_tmp, sr_M_names, sr_M_values_tmp,
    S_stiffness, x, obj_coords;
  PrintStartProc(procname);

  # Substitute -1 entries of coords with the corresponding object length
  obj_coords := [seq(`if`(coords[i] = -1, objs[i]["length"], coords[i]), i = 1..nops(coords))];

  for i from 1 to nops(objs) do
    if IsRod(objs[i]) and
        (ListPadding(obj_coords[i], 3) <> [0, 0, 0]) and
        (ListPadding(eval(obj_coords[i]^~2), 3) <> [eval(objs[i]["length"]^~2), 0, 0]) then
      error "SUPPORT objects can only be applied at extremes of ROD objects"
    end if;
    if IsRod(objs[i]) and (constrained_dof[4..6] <> [0, 0, 0]) then
      error "ROD objects supports can only have translational constraints"
    end if;
  end do;

  if type(stiffness, procedure) then
    S_stiffness := unapply(stiffness(x) *~ constrained_dof, x);
    # Check for non zero stiffness on constrained dof
    if has(map(evalb, stiffness(x)[remove(x -> x=0, ([seq(i, i = 1..6)]) *~ constrained_dof)] <>~ 0), false) then
      error "stiffness corresponding to constrained degrees of freedom cannot be zero";
    end if;
    # Check for zero stiffness on unconstrained dof
    if (S_stiffness(x) <> stiffness(x)) and (not suppress_warnings) then
      WARNING("stiffness components not corresponding to constrained_dof are ignored");
    end if;
  else
    # Check for non zero stiffness on constrained dof
    if has(stiffness[remove(x -> x=0, ([seq(i, i = 1..6)]) *~ constrained_dof)], 0) then
      error "stiffness corresponding to constrained degrees of freedom cannot be zero";
    end if;
    # Check for zero stiffness on unconstrained dof
    S_stiffness := (x) -> stiffness *~ constrained_dof;
    if (S_stiffness(x) <> stiffness) and (not suppress_warnings) then
      WARNING("stiffness components not corresponding to constrained_dof are ignored");
    end if;
  end if;

  S := table({
    "type"                = SUPPORT,
    "constrained_dof"     = constrained_dof,
    "admissible_loads"    = constrained_dof,
    "coordinates"         = [[0,0,0], op(map(ListPadding, obj_coords, 3))],
    "name"                = name,
    "frame"               = RF,
    "targets"             = [earth["name"]] union GetNames(objs),
    "variables"           = [],
    "forces"              = [],
    "moments"             = [],
    "constraint_loads"    = [],
    "support_reactions"   = [], # Expressed in support reference frame
    "stiffness"           = S_stiffness,
    "displacements"       = [],
    "frame_displacements" = LinearAlgebra[IdentityMatrix](4,4)
    });

  # Build the temporary joint
  J_tmp := MakeJoint(name, constrained_dof, [earth, op(objs)], S["coordinates"], RF);

  S["variables"]        := J_tmp["variables"];
  S["forces"]           := J_tmp["forces"];
  S["moments"]          := J_tmp["moments"];
  S["constraint_loads"] := J_tmp["constraint_loads"];

  # Retrieve support force reactions
  sr_F_names := [FX, FY, FZ];
  for i from 1 to nops(S["forces"]) do
    if (S["forces"][i]["target"] = earth["name"]) then
      # Project forces in the support reference frame
      sr_F_values_tmp := Project(S["forces"][i]["components"], ground, S["frame"]);
      for j from 1 to 3 do
        if (sr_F_values_tmp[j] <> 0) then
          S["support_reactions"] := [
            op(S["support_reactions"]),
            sr_F_names[j] = -sr_F_values_tmp[j]
            ];
        end if;
      end do;
      break;
    end if;
  end do;

  # Retrieve support moments reactions
  sr_M_names := [MX, MY, MZ];
  for i from 1 to nops(S["moments"]) do
    if (S["moments"][i]["target"] = earth["name"]) then
      # Project moments in the support reference frame
      sr_M_values_tmp := Project(S["moments"][i]["components"], ground, S["frame"]);
      for j from 1 to 3 do
        if (sr_M_values_tmp[j] <> 0) then
          S["support_reactions"] := [
            op(S["support_reactions"]),
            sr_M_names[j] = -sr_M_values_tmp[j]
            ];
        end if;
      end do;
      break;
    end if;
  end do;

  PrintEndProc(procname);
  return op(S);
end proc: # MakeSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsSupport := proc(
  obj::{anything}, # Object to be checked
  $)::{boolean};

  description "Check if the object <obj> is a SUPPORT object.";

  PrintStartProc(procname);
  PrintEndProc(procname);
  return evalb(obj["type"] = SUPPORT);
end proc: # IsSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsCompliantSupport := proc(
  obj::{anything}, # Object to be checked
  $)::{boolean};

  description "Check if the object <obj> is a SUPPORT object with compliant "
    "constraints.";

  local out, i;
  PrintStartProc(procname);

  out := false;
  if IsSupport(obj) then
    for i from 1 to 6 do
      if (obj["stiffness"](x)[i] <> infinity) and
         (obj["constrained_dof"][i] = 1) then
        out := true;
        break;
      end if;
    end do;
  end if;

  PrintEndProc(procname);
  return out;
end proc; # IsCompliantSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanSupport := proc(
  obj::{SUPPORT}, # Support to be cleaned
  $)::{nothing};

  description "Clean SUPPORT object <obj> internal variables";

  PrintStartProc(procname);
  obj["constraint_loads"]         := [];
  obj["support_reactions"]        := [];
  obj["displacements"]            := [];
  PrintEndProc(procname);
end proc: # CleanSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeJoint := proc(
  name::{string},                                               # Joint name
  constrained_dof::{list},                                      # Constrained degree of freedom
  objs::{list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH})}, # Target objects
  coords::{list},                                               # Joint locations
  RF::{FRAME} := ground,                                        # Reference frame
  {
    stiffness::{procedure,list(algebraic)} := [ # Stiffness components (default = infinite)
      infinity, infinity, infinity,
      infinity, infinity, infinity
    ] *~ constrained_dof,
    shell_objs::{ # Objects to be considered connected to the shell of the joint
      list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH})
    } := [objs[1]]
  }, $)::JOINT;

  description "Make a JOINT object with inputs: joint name <name>, constrained "
    "degrees of freedom <constrained_dof>, target objects <objs>, joint locations "
    "<coords>, and optional reference frame <RF> in which the joint is defined "
    "(default = ground). The optional input <stiffness> is a list of stiffness "
    "components (default = infinite) in the order: [ktx, kty, ktz, krx, kry, krz] "
    " and <shell_objs> is a list of objects to be considered connected to the "
    "shell of the joint (default = [objs[1]]).";

  local J, i, jf_comp, jm_comp, jf_comp_obj, jm_comp_obj, jm_surv, jf_surv,
    jf_comp_cons, jm_comp_cons, constraint, P_tmp, obj_coords, J_stiffness;
  PrintStartProc(procname);

  # Substitute -1 entries of coords with the corresponding object length
  obj_coords := [seq(`if`(coords[i] = -1, objs[i]["length"], coords[i]), i = 1..nops(coords))];

  for i from 1 to nops(objs) do
    if IsRod(objs[i]) then
      # x coordinate of the joint location
      Simplify(ListPadding(eval(obj_coords[i]), 3))[1];
      # x coordinate of the joint location minus target length
      Simplify(eval(%^2) - eval(objs[i]["length"]^~2))[1];
      if  %% <> 0 and %% <> 0. and
          % <> 0 and % <> 0. then
        error "JOINT objects can only be applied at extremes of ROD objects";
      end if;
    end if;
    if IsRod(objs[i]) and (constrained_dof[4..6] <> [0, 0, 0]) then
      error "ROD objects supports can only have translational constraints";
    end if;
  end do;

  if type(stiffness, procedure) then
    J_stiffness := unapply(stiffness(x) *~ constrained_dof, x);
    # Check for non zero stiffness on constrained dof
    if has(map(evalb, stiffness(x)[remove(x -> x=0, ([seq(i, i = 1..6)]) *~ constrained_dof)] <>~ 0), false) then
      error "stiffness corresponding to constrained degrees of freedom cannot be zero";
    end if;
    # Check for zero stiffness on unconstrained dof
    if (J_stiffness(x) <> stiffness(x)) and (not suppress_warnings) then
      WARNING("stiffness components not corresponding to constrained_dof are ignored");
    end if;
  else
    # Check for non zero stiffness on constrained dof
    if has(stiffness[remove(x -> x=0, ([seq(i, i = 1..6)]) *~ constrained_dof)], 0) then
      error "stiffness corresponding to constrained degrees of freedom cannot be zero";
    end if;
    # Check for zero stiffness on unconstrained dof
    J_stiffness := (x) -> stiffness *~ constrained_dof;
    if (J_stiffness(x) <> stiffness) and (not suppress_warnings) then
      WARNING("stiffness components not corresponding to constrained_dof are ignored");
    end if;
  end if;

  J := table({
    "type"                = JOINT,
    "constrained_dof"     = constrained_dof,
    "admissible_loads"    = constrained_dof,
    "coordinates"         = map(ListPadding, obj_coords, 3),
    "name"                = name,
    "frame"               = RF,
    "targets"             = GetNames(objs),
    "shell_targets"       = GetNames(shell_objs),
    "variables"           = [],
    "forces"              = [],
    "moments"             = [],
    "constraint_loads"    = [],
    "stiffness"           = J_stiffness,
    "displacements"       = [],
    "frame_displacements" = LinearAlgebra[IdentityMatrix](4,4)
    });

  # Check if joint position on each object is the same
  # FIXME: does not work for mixed numerical and symbolic coordinates
  #if nops(ells) > 1 then
  #  P_tmp := objs[1]["frame"].Translate(ells[1], 0, 0);
  #  for i from 2 to nops(ells) do
  #    if not (norm(P_tmp - objs[i]["frame"].Translate(ells[i], 0, 0)) = 0) then
  #      error "Joint locations are not the same on all objects";
  #    end if;
  #  end do;
  #end if;

  # Add all the bodies forces
  for i from 1 to nops(objs) do
    # Create joint forces force
    jf_comp := <
      JFx_||(J["name"])||_||(objs[i]["name"]),
      JFy_||(J["name"])||_||(objs[i]["name"]),
      JFz_||(J["name"])||_||(objs[i]["name"])
      >;
    # Keep components compatible with the joint constrained dof
    jf_comp_cons := convert( jf_comp *~ <op(constrained_dof[1..3])>,
      list);
    # Extract the survived components
    jf_surv := remove(x -> x = 0, jf_comp_cons);
    # Project the components into object frame and extract admissible loads
    jf_comp_obj := convert(
      Project(jf_comp_cons, RF, objs[i]["frame"])
      .~ <op(objs[i]["admissible_loads"][1..3])>,
      list);
    # Check if there are reactions
    if (nops(jf_surv) <> 0) then
      # Create the reaction force between joint and obj
      # Force on obj
      JF_||(name)||_||(objs[i]["name"]) := MakeForce(
        jf_comp_obj, obj_coords[i], objs[i], objs[i]["frame"]
        );
      # Force on joint
      JF_||(objs[i]["name"])||_||(name) := MakeForce(
        -jf_comp_cons, 0, J, RF);
      # Use the non admissible loads to build the loads constraint
      constraint := convert(
        Project(jf_comp_cons, RF, objs[i]["frame"])
        .~ <eval(map((x->evalb(x = 0)), objs[i]["admissible_loads"][1..3]), [true = 1, false = 0])>,
        list);
      # Remove the null equations
      constraint := remove(x -> x = 0, constraint);
      # Update the joint constraint loads
      J["constraint_loads"] := J["constraint_loads"] union constraint;
      # Update the output joint
      J["variables"] := J["variables"] union jf_surv;
      J["forces"] := J["forces"] union
       [JF_||(name)||_||(objs[i]["name"]),
        JF_||(objs[i]["name"])||_||(name)];
    end if;
  end do;

  # Add all the bodies moments
  for i from 1 to nops(objs) do
    # Create moment compatible with joint constrained dof
    jm_comp := <
      JMx_||(J["name"])||_||(objs[i]["name"]),
      JMy_||(J["name"])||_||(objs[i]["name"]),
      JMz_||(J["name"])||_||(objs[i]["name"])
      >;
    # Keep components compatible with the joint constrained dof
    jm_comp_cons := convert(jm_comp *~ <op(constrained_dof[4..6])>,
      list);
    # Extract the survived components
    jm_surv := remove(x -> x = 0, jm_comp_cons);
    # Project the components into object frame and extract the admissible loads
    jm_comp_obj := convert(
      Project(jm_comp_cons, RF, objs[i]["frame"])
      .~ <op(objs[i]["admissible_loads"][4..6])>,
      list);
    # Check if there are reactions
    if (nops(jm_surv) <> 0) then
      # Create the reaction force between joint and obj
      # Moment on obj
      JM_||(name)||_||(objs[i]["name"]) := MakeMoment(
        jm_comp_obj, obj_coords[i], objs[i], objs[i]["frame"]
        );
      # Moment on joint
      JM_||(objs[i]["name"])||_||(name) := MakeMoment(
      -jm_comp_cons, 0, J, RF);
      # Use the non admissible loads to build the loads constraint
      constraint := convert(
        Project(jm_comp_cons, RF, objs[i]["frame"])
        .~ <eval(map((x->evalb(x = 0)), objs[i]["admissible_loads"][4..6]), [true = 1, false = 0])>,
        list);
      constraint := remove(x -> x = 0, constraint);
      # Update the joint constraint loads
      J["constraint_loads"] := J["constraint_loads"] union constraint;
      # Update the output joint
      J["variables"] := J["variables"] union jm_surv;
      J["moments"] := J["moments"] union
        [JM_||(name)||_||(objs[i]["name"]),
        JM_||(objs[i]["name"])||_||(name)];
    end if;
  end do;

  PrintEndProc(procname);
  return op(J);
end proc: # MakeJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsJoint := proc(
  obj::{anything}, # Object to be checked
  $)::{boolean};

  description "Check if the objecy <obj> is a JOINT object.";

  PrintStartProc(procname);
  PrintEndProc(procname);
  return evalb(obj["type"] = JOINT);
end proc: # IsJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsCompliantJoint := proc(
  obj::{anything}, # Object to be checked
  $)::{boolean};

  description "Check if the object <obj> is a JOINT object with compliant "
    "constraints.";

  local out, i;
  PrintStartProc(procname);

  out := false;
  if IsJoint(obj) then
    for i from 1 to 6 do
      if (obj["stiffness"](x)[i] <> infinity) and
         (obj["constrained_dof"][i] = 1) then
        out := true;
        break;
      end if;
    end do;
  end if;

  PrintEndProc(procname);
  return out;
end proc; # IsCompliantJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanJoint := proc(
  obj::{JOINT}, # Object to be cleaned
  $)::{nothing};

  description "Clean JOINT object <obj> internal variables.";

  PrintStartProc(procname);
  obj["constraint_loads"] := [];
  PrintEndProc(procname);
end proc: # CleanJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeRodPoints := proc(
  name::{string},  # Object name
  point1::{POINT}, # First point
  point2::{POINT}, # Second point
  vec::{VECTOR},   # Vector for XY-plane
  {
    area::{algebraic, procedure} := infinity,      # Section area (m^2)
    material::{MATERIAL}         := MakeMaterial() # Material
  }, $)::{ROD};

  description "Create a ROD object with inputs: object name <name>, first "
    "point <point1>, second point <point2>, vector in XY-plane <vec>, optional "
    "section area <area> and material type <material>.";

  local ell, ex, ey, ez, RF, out;
  PrintStartProc(procname);

  if (point1 = point2) then
    error "input points are the same";
  end if;

  # FIXME: does not work with symbolic vectors
  # if (Norm2(vec) < 0) then
  #   error "input vector is null";
  # end if;

  ell := Norm2(point2 - point1);
  ex  := (point2 - point1) /~ ell;
  ey  := convert(LinearAlgebra[CrossProduct](<op(vec)>, <op(ex)>), list);
  ey  := ey /~ Norm2(ey);
  ez  := convert(LinearAlgebra[CrossProduct](<op(ex)>, <op(ey)>), list);
  ez  := ez /~ Norm2(ez);

  RF := <<ex[1],     ex[2],     ex[3],     0>|
         <ey[1],     ey[2],     ey[3],     0>|
         <ez[1],     ez[2],     ez[3],     0>|
         <point1[1], point1[2], point1[3], 1>>;

  out := MakeRod(
    name, ell, Simplify(RF),
    parse("area")     = area,
    parse("material") = material
    );

  PrintEndProc(procname);
  return op(out);
end proc: # MakeRodPoints

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeRod := proc(
  name::{string},        # Object name
  ell::{algebraic},      # Length (m)
  RF::{FRAME} := ground, # Reference frame
  {
    area::{algebraic, procedure} := infinity,      # Section area (m^2)
    material::{MATERIAL}         := MakeMaterial() # Material
  }, $)::{ROD};

  description "Create a ROD object with inputs: object name <name>, reference "
    "length <ell>, optional reference frame <RF> in which the rod is defined, "
    "and optional section area <area> and material type <material>.";

  local area_proc, out;
  PrintStartProc(procname);

  if type(area, procedure) then
    area_proc := area;
  else
    area_proc := (x) -> area;
  end if;

  out := table({
    "type"                = ROD,
    "name"                = name,
    "length"              = ell,
    "area"                = area_proc,
    "material"            = material,
    "frame"               = RF,
    "admissible_loads"    = [1, 0, 0, 0, 0, 0],
    "internal_actions"    = [],
    "displacements"       = [],
    "frame_displacements" = LinearAlgebra[IdentityMatrix](4,4)
    });

  PrintEndProc(procname);
  return op(out);
end proc: # MakeRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsRod := proc(
  obj::{anything}, # Object to be checked
  $)::{boolean};

  description "Check if the object <obj> is a ROD object.";

  PrintStartProc(procname);
  PrintEndProc(procname);
  return evalb(obj["type"] = ROD);
end proc: # IsRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanRod := proc(
  obj::{ROD}, # Object to be cleaned
  $)::{nothing};

  description "Clean ROD object <obj> internal variables.";

  PrintStartProc(procname);
  obj["internal_actions"] := [];
  obj["displacements"]    := [];
  PrintEndProc(procname);
end proc: # CleanRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeBeamPoints := proc(
  name::{string},  # Object name
  point1::{POINT}, # First point
  point2::{POINT}, # Second point
  vec::{VECTOR},   # Vector for XY-plane (normal vector)
  {
    area::{algebraic, procedure}                   := infinity,       # Section area (m^2)
    timo_shear_coeff::{list(algebraic), procedure} := [5/6, 5/6],     # Timoshenko shear coefficient
    material::{MATERIAL}                           := MakeMaterial(), # Material object
    I_xx::{algebraic, procedure}                   := infinity,       # Section x-axis inertia (m^4)
    I_yy::{algebraic, procedure}                   := infinity,       # Section y-axis inertia (m^4)
    I_zz::{algebraic, procedure}                   := infinity        # Section z-axis inertia (m^4)
  }, $)::{BEAM};

  description "Create a BEAM object with inputs: object name <name>, first "
    "point <point1>, second point <point2>, vector in XY-plane <vec>,  optional "
    "section area <area>, optional Timoshenko shear coefficient <timo_shear_coeff>, "
    "optional material type <material>, optional section x-axis inertia <I_xx>, "
    "optional section y-axis inertia <I_yy> and optional section z-axis inertia "
    "<I_zz>.";

  local ell, ex, ey, ez, RF, out;
  PrintStartProc(procname);

  if (point1 = point2) then
    error "input points are the same";
  end if;

  # FIXME: does not work with symbolic vectors
  # if (Norm2(vec) < 0) then
  #   error "input vector is null";
  # end if;

  ell := Norm2(point2 - point1);
  ex  := (point2 - point1) /~ ell;
  ey  := convert(LinearAlgebra[CrossProduct](<op(vec)>, <op(ex)>), list);
  ey  := ey /~ Norm2(ey);
  ez  := convert(LinearAlgebra[CrossProduct](<op(ex)>, <op(ey)>), list);
  ez  := ez /~ Norm2(ez);


  RF := <<ex[1],     ex[2],     ex[3],     0>|
         <ey[1],     ey[2],     ey[3],     0>|
         <ez[1],     ez[2],     ez[3],     0>|
         <point1[1], point1[2], point1[3], 1>>;

  out := MakeBeam(name, ell, Simplify(RF),
    parse("area")             = area,
    parse("timo_shear_coeff") = timo_shear_coeff,
    parse("material")         = material,
    parse("I_xx")             = I_xx,
    parse("I_yy")             = I_yy,
    parse("I_zz")             = I_zz
    );

  PrintEndProc(procname);
  return op(out);
end proc: # MakeBeamPoints

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeBeam := proc(
  name::{string},        # Object name
  ell::{algebraic},      # Length (m)
  RF::{FRAME} := ground, # Reference frame
  {
    area::{algebraic, procedure}                   := infinity,       # Section area (m^2)
    timo_shear_coeff::{list(algebraic), procedure} := [5/6, 5/6],     # Timoshenko shear coefficient
    material::{MATERIAL}                           := MakeMaterial(), # Material object
    I_xx::{algebraic, procedure}                   := infinity,       # Section x-axis inertia (m^4)
    I_yy::{algebraic, procedure}                   := infinity,       # Section y-axis inertia (m^4)
    I_zz::{algebraic, procedure}                   := infinity        # Section z-axis inertia (m^4)
  }, $)::{BEAM};

  description "Create a BEAM object with inputs: object name <name>, reference "
    "length <ell>, optional reference frame <RF> in which the rod is defined, and "
    "optional section area <area> and material type <material> and inertias on "
    "x- <I_xx>, y- <I_yy>, and z-axis <I_zz>.";

  local area_proc, timo_shear_coeff_proc, I_xx_proc, I_yy_proc, I_zz_proc, out;
  PrintStartProc(procname);

  if type(area, procedure) then
    area_proc := area;
  else
    area_proc := (x) -> area;
  end if;

  if type(timo_shear_coeff, procedure) then
    timo_shear_coeff_proc := timo_shear_coeff;
  else
    timo_shear_coeff_proc := (x) -> timo_shear_coeff;
  end if;

  if type(I_xx, procedure) then
    I_xx_proc := I_xx;
  else
    I_xx_proc := (x) -> I_xx;
  end if;

  if type(I_yy, procedure) then
    I_yy_proc := I_yy;
  else
    I_yy_proc := (x) -> I_yy;
  end if;

  if type(I_zz, procedure) then
    I_zz_proc := I_zz;
  else
    I_zz_proc := (x) -> I_zz;
  end if;

  out := table({
    "type"                = BEAM,
    "name"                = name,
    "length"              = ell,
    "area"                = area_proc,
    "timo_shear_coeff"    = timo_shear_coeff_proc,
    "material"            = material,
    "inertias"            = [I_xx_proc, I_yy_proc, I_zz_proc],
    "frame"               = RF,
    "admissible_loads"    = [1, 1, 1, 1, 1, 1],
    "internal_actions"    = [],
    "displacements"       = [],
    "frame_displacements" = LinearAlgebra[IdentityMatrix](4,4)
    });

  PrintEndProc(procname);
  return op(out);
end proc: # MakeBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsBeam := proc(
  obj::{anything}, # Object to be checked
  $)::{boolean};

  description "Check if the object <obj> is a BEAM object.";

  PrintStartProc(procname);
  PrintEndProc(procname);
  return evalb(obj["type"] = BEAM);
end proc: # IsBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanBeam := proc(
  obj::{BEAM}, # Object to be cleaned
  $)::{nothing};

  description "Clean BEAM object <obj> internal variables.";

  PrintStartProc(procname);
  obj["internal_actions"] := [];
  obj["displacements"]    := [];
  PrintEndProc(procname);
end proc: # CleanBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeRigidBody := proc(
  name::{string},        # Object name
  RF::{FRAME} := ground, # Reference frame
  {
    COM::{list(algebraic)} := [0, 0, 0], # COM position in RF (default: [0, 0, 0])
    mass::{algebraic}      := 0          # Mass (kg)
  },
  $)::{RIGID_BODY};

  description "Create a RIGID_BODY object with inputs: object name <name>, "
    "reference frame <RF> in which the rigid body is defined, and optional "
    "center of mass position <COM> and mass <mass>.";

  local out;
  PrintStartProc(procname);

  out := table({
    "type"                = RIGID_BODY,
    "name"                = name,
    "frame"               = RF,
    "COM"                 = COM,
    "mass"                = mass,
    "admissible_loads"    = [1, 1, 1, 1, 1, 1],
    "frame_displacements" = LinearAlgebra[IdentityMatrix](4,4)
    });

  PrintEndProc(procname);
  return op(out);
end proc: # MakeRigidBody

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsRigidBody := proc(
  obj::{anything}, # Object to be checked
  $)::{boolean};

  description "Check if the object <obj> is a RIGID_BODY object.";

  PrintStartProc(procname);
  PrintEndProc(procname);
  return evalb(obj["type"] = RIGID_BODY);
end proc: # IsRigidBody

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeSpringDisplacement := proc(
  spring_load::{algebraic},      # Load on the spring
  spring_stiffness::{procedure}, # Spring stiffness
  $)::{algebraic};

  description "Compute the displacement of a spring give the load <spring_load> "
  "and spring stiffness <stiffness>.";

  local x, out, Dx;
  PrintStartProc(procname);

  #Physics[Assume](spring_load * Dx > 0);

  out := RealDomain[solve](
    spring_load = convert(integrate(spring_stiffness(x), x = 0..Dx), signum), Dx
    );

  if nops([out]) > 1 then
    out := convert([out], piecewise);
  else
    out := convert(out, piecewise);
  end if;

  PrintEndProc(procname);
  return out;
end proc: # ComputeSpringDisplacement

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeSpringEnergy := proc(
  spring_load::{algebraic},      # Load on the spring
  spring_stiffness::{procedure}, # Spring stiffness
  $)::{algebraic};

  description "Compute the potential energy of a spring give the load <spring_load> "
  "and spring stiffness <stiffness>.";

  local disp, x, out;
  PrintStartProc(procname);

  disp := ComputeSpringDisplacement(spring_load, spring_stiffness);
  out  := integrate(integrate(spring_stiffness(x), x), x = 0..disp);

  PrintEndProc(procname);
  return out;
end proc: # ComputeSpringEnergy

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeSupportDisplacements := proc(
  obj::{SUPPORT}, # Support object
  $)::{nothing};

  description "Compute the displacements of the support <obj> from its "
    "support reactions.";

  local sup_disp, disp_vec, i, disp, x, sup_reac;
  PrintStartProc(procname);

  sup_disp := [];
  disp_vec := [tx, ty, tz, rx, ry, rz];
  sup_reac := ['FX', 'FY', 'FZ', 'MX', 'MY', 'MZ'];

  for i from 1 to 6 do
    if (obj["constrained_dof"][i] = 1) and
        member(sup_reac[i], map(lhs, obj["support_reactions"])) then
      disp := ComputeSpringDisplacement(subs(obj["support_reactions"], - sup_reac[i]),
        (x -> obj["stiffness"](x)[i]));
      sup_disp := sup_disp union [disp_vec[i] = disp];
    end if;
  end do;

  if (verbose_mode > 0) then
    printf(
      "%*sMessage (in ComputeSupportDisplacements) updating %s %s's displacements...\n",
      print_indent, "|   ", obj["type"], obj["name"]
      );
  end if;

  obj["displacements"] := sup_disp;

  if (verbose_mode > 0) then
    printf("%*sDONE\n", print_indent, "|   ");
  end if;

  PrintEndProc(procname);
end proc: # ComputeSupportDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeJointDisplacements := proc(
  obj::{JOINT},     # Joint object
  sol::{list, set}, # List of solutions for joint forces
  $)::{nothing};

  description "Compute the displacements of the joint <obj>.";

  local jnt_disp, disp_vec, i, disp, x, jnt_load, f;
  PrintStartProc(procname);

  jnt_disp := [];
  disp_vec := [tx, ty, tz, rx, ry, rz];
  jnt_load := [0,0,0,0,0,0];

  for i from 1 to 6 do
    if (obj["constrained_dof"][i] = 1) then
        if i<4 then
          # Forces
          for f in obj["forces"] do
            if member(f["target"], obj["shell_targets"]) then
              jnt_load[i] := jnt_load[i] + f["components"][i];
            end if;
          end do;
        else
          # Moments
          for f in obj["moments"] do
            if member(f["target"], obj["shell_targets"]) then
              jnt_load[i] := jnt_load[i] + f["components"][i-3];
            end if;
          end do;
        end if;
      disp := ComputeSpringDisplacement(Subs(sol, jnt_load[i]),
        (x -> obj["stiffness"](x)[i]));
      jnt_disp := jnt_disp union [disp_vec[i] = disp];
    end if;
  end do;

  if (verbose_mode > 0) then
    printf(
      "%*sMessage (in ComputeJointDisplacements) updating %s %s's displacements...\n",
      print_indent, "|   ", obj["type"], obj["name"]
      );
  end if;

  obj["displacements"] := jnt_disp;

  if (verbose_mode > 0) then
    printf("%*sDONE\n", print_indent, "|   ");
  end if;

  PrintEndProc(procname);
end proc: # ComputeJointDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeStructure := proc(
  objs::{ # Structure objects
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set( {BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  exts::{ # External actions
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  } := [],
  {
    hyper_vars::{list ,set} := [],
      # Hyperstatic variables
    hyper_disp::{list, set} := [seq(0, 1..nops(hyper_vars))]
      # Hyperstatic displacements
  }, $)::{STRUCTURE};

  description "Create a STRUCTURE object with inputs: structure objects <objs>, "
    "external actions <exts>, optional hyperstatic variables <hyper_vars>, "
    "optional hyperstatic displacements <hyper_disp>.";

  local num_dof, names, candidate_hyp_vars, Graph, obj, S_ext, out;
  PrintStartProc(procname);

  # Check for duplicate names
  names := [];
  for obj in objs do
    if member(obj["name"], names) then
      error "duplicate names found on structure objects";
    end if;
    names := names union [obj["name"]];
  end do;

  num_dof, Graph := ComputeDOF(objs);

  if (num_dof < 0) then
    if (nops(hyper_vars) <> -num_dof) and (not suppress_warnings) then
      candidate_hyp_vars := [];
      for obj in objs do
        if IsSupport(obj) or IsJoint(obj) then
          candidate_hyp_vars := candidate_hyp_vars union obj["variables"];
        end if;
      end do;
    WARNING(
      "the structure is hyperstatic with %1 overconstrained directions, "
      "please check the structure supports and joints. Also consider defining "
      "the hyperstatic variables by adding 'hyper_vars' property in the "
      "'MakeStructure' method or simply defining 'hyperstatic_variables' "
      "field in an already existing STRUCTURE object by choosing from the "
      "folloving hyperstatic candidate variables: %2.",
      abs(num_dof), candidate_hyp_vars
      );
    else
      if (verbose_mode > 0) then
        printf("%*sMessage (in MakeStructure) "
          "hyperstatic structure detected with %d overconstrained directions\n",
          print_indent, "|   ", abs(num_dof));
      end if;
    end if;
  elif (num_dof > 0) and (not suppress_warnings) then
    #error "not enough constraints in the structure";
    WARNING(
      "the structure is underconstrained with %1 unconstrained directions. "
      "Results computation may fail due to rigid body motions.",
      num_dof
    );
  else
    if (verbose_mode > 0) then
      printf("%*sMessage (in MakeStructure) isostatic structure detected\n", print_indent, "|   ");
    end if;
  end if;

  # Add gravity distributed load
  S_ext := exts;
  if (gravity <> [0, 0, 0]) then
    for obj in objs do
      if IsRod(obj) and (not suppress_warnings) then
        WARNING("Message (in SolveStructure) gravity load is not supported for rod %1", obj);
      elif IsBeam(obj) then
        g_load||(obj["name"]) := MakeQForce(
          (x -> gravity *~ obj["area"](x) *~ obj["material"]["density"]),
          obj,ground
          );
        S_ext := S_ext union {g_load||(obj["name"])};
      elif IsRigidBody(obj) then
        g_load||(obj["name"]) := MakeForce(
          gravity *~ obj["mass"], obj["COM"], obj, ground
          );
        S_ext := S_ext union {g_load||(obj["name"])};
      end if;
    end do;
  end if;

  out := table({
    "type"                       = STRUCTURE,
    "objects"                    = objs,
    "external_actions"           = S_ext,
    "dof"                        = num_dof,
    "connections_graph"          = Graph,
    "hyperstatic_variables"      = convert(hyper_vars, list),
    "hyperstatic_displacements"  = convert(hyper_disp, list),
    "equations"                  = [],
    "variables"                  = [],
    "potential_energy"           = NULL,
    "veils"                      = [],
    "support_reactions_solved"   = false,
    "internal_actions_solved"    = false,
    "displacements_solved"       = false,
    "potential_energy_solved"    = false,
    "frame_displacements_solved" = false
    });

  PrintEndProc(procname);
  return out;
end proc: # MakeStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsStructure := proc(
  obj::{anything}, # Object to be checked
  $)::{boolean};

  description "Check if the object <obj> is a STRUCTURE object.";

  PrintStartProc(procname);
  PrintEndProc(procname);
  return evalb(obj["type"] = STRUCTURE);
end proc: # IsStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanStructure := proc(
  obj::{STRUCTURE}, # Object to be cleaned
  $)::{nothing};

  description "Clean STRUCTURE object <obj> internal variables.";

  local i;
  PrintStartProc(procname);

  # Clean internal variables
  obj["equations"]                  := [];
  obj["variables"]                  := [];
  obj["potential_energy"]           := NULL;
  obj["veils"]                      := [];
  obj["support_reactions_solved"]   := false;
  obj["internal_actions_solved"]    := false;
  obj["displacements_solved"]       := false;
  obj["potential_energy_solved"]    := false;
  obj["frame_displacements_solved"] := false;

  # Clean objects
  for i from 1 to nops(obj["objects"]) do
    if IsBeam(obj[i]) then
      obj["objects"][i] := CleanBeam(i);
    elif IsRod(obj[i]) then
      obj["objects"][i] := CleanRod(i);
    elif IsSupport(obj[i]) then
      obj["objects"][i] := CleanSupport(i);
    elif IsJoint(obj[i]) then
      obj["objects"][i] := CleanJoint(i);
    end if;
  end do;

  PrintEndProc(procname);
end proc: # CleanStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeDOF := proc(
  objs::{ # Structure objects
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set( {BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  }, $)::integer ::function;

  description "Compute the degree of freedom of the input structure objects "
    "<objs>.";

  local dof, objs_tmp, obj, i, j, k, vertex, colors, G;
  PrintStartProc(procname);

  dof      := 0;
  objs_tmp := objs union [earth];

  # Built connections graph
  vertex := [];
  colors := [];
  if (verbose_mode > 0) then
    printf("%*sMessage (in ComputeDOF) checking structure connections...\n", print_indent, "|   ");
  end if;
  for i from 1 to nops(objs_tmp) do
    vertex := vertex union [objs_tmp[i]["name"]];
    colors := colors union [ObjectColor(objs_tmp[i])];
  end do;
  G := GraphTheory[Graph](vertex);
  GraphTheory[HighlightVertex](G, vertex, colors);
  for i from 1 to nops(objs_tmp) do
    if IsSupport(objs_tmp[i]) or IsJoint(objs_tmp[i]) then
      for j from 1 to nops(objs_tmp) do
        if (member(objs_tmp[j]["name"], objs_tmp[i]["targets"])) then
          GraphTheory[AddEdge](G, {objs_tmp[i]["name"], objs_tmp[j]["name"]});
        end if;
      end do;
    end if;
  end do;

  if (verbose_mode > 0) then
    printf("%*sDONE\n", print_indent, "|   ");
    printf("%*sMessage (in ComputeDOF) display connections graph...\n", print_indent, "|   ");
    if (verbose_mode > 1) then
      print(GraphTheory[DrawGraph](G), layout = tree);
    end if;
  end if;

  # Check graph connections
  if GraphTheory[IsConnected](G) then
    if (verbose_mode > 0) then
      printf("%*sDONE\n", print_indent, "|   ");
    end if;
  else
    print(GraphTheory[DrawGraph](G), layout = tree);
    WARNING("unconnected elements detected in the structure");
  end if;

  if (verbose_mode > 0) then
    printf("%*sMessage (in ComputeDOF) computing degrees of freedom...\n", print_indent, "|   ");
  end if;

  for obj in objs_tmp do
    if IsBeam(obj) then
      dof := dof + 6;
    elif IsRigidBody(obj) then
      dof := dof + 6;
    elif IsRod(obj) then
      dof := dof + 5;
    elif IsJoint(obj) then
      dof := dof - add(obj["constrained_dof"][k], k = 1..6) * (nops(obj["targets"]) - 1);
    elif IsSupport(obj) then
      dof := dof - add(obj["constrained_dof"][k], k = 1..6) * (nops(obj["targets"]) - 1);
    end if;
  end do;

  if (verbose_mode > 0) then
    printf("%*sDONE - DOF = %d\n", print_indent, "|   ", dof);
  end if;

  PrintEndProc(procname);
  if _nresults = 2 then
    return dof, G;
  else
    return dof;
  end
end proc: # ComputeDOF

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

DrawStructureGraph := proc(
  obj::{STRUCTURE}, # Object to be cleaned
  $)::{procedure};

  description "Draw the connections graph of the STRUCTURE object <obj>.";

  local  out;
  PrintStartProc(procname);

  out := plots[display](
    GraphTheory[DrawGraph](obj["connections_graph"], layout = tree),
    title = "Structure connections graph");

  PrintEndProc(procname);
  return out;
end proc: # DrawStructureGraph

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

DrawStructureSparseMatrix := proc(
  obj::{STRUCTURE}, # Object to be cleaned
  {
    GaussianElimination::{boolean} := false # Apply Gaussian elimination
  },
  $)::{procedure};

  description "Draw the sparse matrix for the equation system of STRUCTURE "
    "object <obj>.";

  local  out, A, b;
  PrintStartProc(procname);

  A,b := LinearAlgebra[GenerateMatrix](obj["equations"], obj["variables"]);
  if (GaussianElimination) then
    A := LinearAlgebra["GaussianElimination"](A);
  end if;
  out := plots[display](
    plots[sparsematrixplot](A,matrixview),
    title = "Structure sparse matrix");

  PrintEndProc(procname);
  return out;
end proc: # DrawStructureSparseMatrix

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

NewtonEuler := proc(
  exts::{ # External actions
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  },
  obj::{BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}, # Object to compute the equilibrium
  {
    pole::{POINT} := [0, 0, 0],             # Pole to compute the equilibrium
    upper_lim::{algebraic} := obj["length"] # Upper limit of the integration
  }, $)::{list};

  description "Compute the Newton-Euler static equilibrium equations given a set "
    "of external actions <exts>, and object to compute the equilibrium <obj>, "
    "the axial coordinate of the pole <pole>, and an optional upper limit of the "
    "integration <upper_lim>.";

  local eq_T, eq_R, i, x, arm, out;
  PrintStartProc(procname);

  eq_T := [0, 0, 0];
  for i from 1 to nops(exts) do
    if exts[i]["target"] = obj["name"] then
      if IsForce(exts[i]) then
        eq_T := eq_T + exts[i]["components"];
      elif IsQForce(exts[i]) then
        eq_T := eq_T + map(
          integrate, exts[i]["components"](x), x = 0..upper_lim
          );
      end if;
    elif (not suppress_warnings) then
      WARNING("Message (in NewtonEuler) %1 is not applied to %2", exts[i], obj);
    end if;
  end do;

  eq_R := [0, 0, 0];
  for i from 1 to nops(exts) do
    if (exts[i]["target"] = obj["name"]) then
      if IsMoment(exts[i]) then
        eq_R := eq_R + exts[i]["components"];
      elif IsForce(exts[i]) then
          arm := <op(pole - exts[i]["coordinate"])>;
        eq_R := eq_R +
          convert(LinearAlgebra[CrossProduct](<op(exts[i]["components"])>, arm), list);
      elif IsQForce(exts[i]) then
        arm := <op(pole - [x, 0, 0])>;
        eq_R := eq_R + map(integrate,
          convert(LinearAlgebra[CrossProduct](<op(exts[i]["components"](x))>, arm), list),
          x = 0..upper_lim);
      elif IsQMoment(FMQ[i]) then
        eq_R := eq_R + map(integrate,
          exts[i]["components"](x), x = 0..upper_lim);
      end if;
    elif (not suppress_warnings) then
      WARNING("Message (in NewtonEuler) %1 is not applied to %2", exts[i], obj);
    end if;
  end do;

  out := eq_T union eq_R;

  PrintEndProc(procname);
  return out;
end proc: # NewtonEuler

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SolveStructure := proc(
  struct::{STRUCTURE}, # Structure object
  {
    compute_internal_actions::{boolean}    := false, # Internal actions computation flag
    compute_displacements::{boolean}       := false, # Displacement computation flag
    compute_potential_energy::{boolean}    := false, # Potential energy computation flag
    compute_frame_displacements::{boolean} := false, # Frame displacements computation flag
    timoshenko_beam::{boolean}             := false, # Timoshenko beam flag
    implicit::{boolean}                    := false, # Implicit solution flag
    unveil_results::{boolean}              := true,  # Unveil results flag
    dummy_vars::{list, set}                := []     # Dummy variables
  }, $)::{STRUCTURE};

  description "Solve the static equilibrium of a structure with inputs: "
    "structure <struct>, optional compute internal action enabling flag "
    "<compute_internal_actions>, optional compute displacement enabling flag "
    "<compute_displacements>, optional Timoshenko beam flag <timoshenko_beam>, "
    "optional implicit solution flag <implicit>, and optional unveil results "
    "flag <unveil_results>.";

  local g_load, S_obj, S_rigid, S_ext, S_support, S_joint, S_con_forces, vars,
    sol, obj, x, str_eq, str_vars, P_energy, veiling_idx, veils;
  PrintStartProc(procname);

  # Clean structure
  CleanStructure(struct);

  # Set veiling_label
  veiling_idx := 1;
  veiling_label := cat('_V', veiling_idx);

  # Parsing inputs
  S_obj        := {};
  S_rigid      := {};
  S_ext        := {};
  S_support    := {};
  S_joint      := {};
  S_con_forces := {};
  vars         := [];
  for obj in struct["objects"] do
    if IsBeam(obj) or IsRod(obj) then
      S_obj := S_obj union {eval(obj)};
    elif IsRigidBody(obj) then
      S_rigid := S_rigid union {eval(obj)};
    elif IsSupport(obj) then
      S_support    := S_support union {eval(obj)};
      S_con_forces := S_con_forces union obj["forces"] union obj["moments"];
      vars         := vars union obj["variables"];
    elif IsJoint(obj) then
      S_joint      := S_joint union {eval(obj)};
      S_con_forces := S_con_forces union obj["forces"] union obj["moments"];
      vars         := vars union obj["variables"];
    end if;
    unassign('obj'); # FIXME: try to remove this line when the bug is fixed
  end do;

  S_ext := struct["external_actions"];

  # Set module local variable keep_veiled
  keep_veiled := not unveil_results;

  # Solve isostatic structure
  if (struct["dof"] >= 0) then

    if struct["dof"] > 0 and (not suppress_warnings) then
      WARNING("Message (in SolveStructure) structure is underconstrained. "
        "Trying to solve it anyway. Results computation may fail due to rigid "
        "body motions.");
    end if;

    if (verbose_mode > 0) then
      printf("%*sMessage (in SolveStructure) solving the isostatic structure...\n", print_indent, "|   ");
    end if;
    sol, str_eq, str_vars := IsostaticSolver(
      S_obj union S_rigid union S_joint union S_support,
      S_ext union S_con_forces,
      vars,
      parse("implicit") = implicit
      );
    # Update Structure equations and variables
    struct["equations"] := str_eq;
    struct["variables"] := str_vars;
    if (verbose_mode > 0) then
      printf("%*sDONE\n", print_indent, "|   ");
      if (verbose_mode > 1) then
        printf("%*sMessage (in SolveStructure) updating support reactions fields...\n", print_indent, "|   ");
      end if;
    end if;
    # Update support reactions properties
    for obj in S_support do
      obj["support_reactions"] := [
        seq(lhs(obj["support_reactions"][j]) = Subs(sol, rhs(obj["support_reactions"][j])),
        j = 1..nops(obj["support_reactions"]))
      ];
    end do;
    if (verbose_mode > 0) then
      printf("%*sDONE\n", print_indent, "|   ");
    end if;

  # Solve hyperstatic structure
  elif (struct["dof"] < 0) then

    if (nops(struct["hyperstatic_variables"]) <> -struct["dof"]) then
      error "mismatch in the structure degrees of freedom, check the hyper"
        "static variables of the structure and update the structure object";
    end if;
    if (verbose_mode > 0) then
      printf("%*sMessage (in SolveStructure) solving the hyperstatic structure\n", print_indent, "|   ");
    end if;
    sol, str_eq, str_vars, P_energy := HyperstaticSolver(
      S_obj union S_rigid union S_joint union S_support,
      S_ext union S_con_forces,
      vars,
      struct["hyperstatic_variables"],
      struct["hyperstatic_displacements"],
      parse("timoshenko_beam") = timoshenko_beam,
      parse("implicit")        = implicit
      );
    if (verbose_mode > 0) then
      printf("%*sDONE\n", print_indent, "|   ");
      if (verbose_mode > 1) then
        printf("%*sMessage (in SolveStructure) hyperstatic solver solution:\n", print_indent, "|   ");
        print(<sol>);
      end if;
    end if;
    # Update Structure equations and variables
    struct["equations"] := str_eq;
    struct["variables"] := str_vars;
    # Update structure energy
    struct["potential_energy"] := P_energy;
    # Set potential energy computed flag
    struct["potential_energy_solved"] := true;
    # Update objects internal actions
    for obj in S_obj do
      obj["internal_actions"] := subs(sol, obj["internal_actions"]);
    end do;
    # Set internal actions computed flag
    struct["internal_actions_solved"] := true;
    # Update support reactions properties
    if (verbose_mode > 0) then
      printf("%*sMessage (in SolveStructure) updating support reactions fields...\n", print_indent, "|   ");
    end if;
    for obj in S_support do
    obj["support_reactions"] := Subs(sol, obj["support_reactions"]);
    end do;
    if (verbose_mode > 0) then
      printf("%*sDONE\n", print_indent, "|   ");
    end if;
  end if;

  # Add veils
  if keep_veiled and type(sol[-1], list) then
    struct["veils"] := struct["veils"] union sol[-1];
    veiling_idx     := veiling_idx + 1;
    veiling_label   := cat('_V', veiling_idx);
  end if;

  # Set support reactions solved flag
  struct["support_reactions_solved"] := true;

  # Compute internal actions
  if ((compute_internal_actions) or (compute_displacements) or (compute_potential_energy)) and not struct["internal_actions_solved"] then
    ComputeInternalActions(
      S_obj, S_ext union S_con_forces, sol
      );
    # Set internal actions computed flag
    struct["internal_actions_solved"] := true;
  end if;

  # Compute potential energy
  if (compute_potential_energy) and not struct["potential_energy_solved"] then
    if implicit then
      error "potential energy cannot be computed in implicit mode";
    end if;
    P_energy := ComputePotentialEnergy(
      S_obj union S_support union S_joint, sol,
      parse("timoshenko_beam") = timoshenko_beam,
      parse("dummy_vars")      = dummy_vars);
    # Update structure energy
    struct["potential_energy"] := P_energy;
    # Set potential energy computed flag
    struct["potential_energy_solved"] := true;
  end if;

  # Compute displacements
  if (compute_displacements) and (not struct["displacements_solved"]) then
    ComputeDisplacements(
      S_obj union S_joint union S_support, S_ext union S_con_forces, sol
      );
    # Set displacements computed flag
    struct["displacements_solved"] := true;
  end if;

  if (compute_frame_displacements) and (not struct["frame_displacements_solved"]) then
    veils := ComputeObjectFrameDisplacements(struct,
                                              parse("timoshenko_beam") = timoshenko_beam,
                                              parse("unveil_results")  = unveil_results);
    # Set frame displacements computed flag
    struct["frame_displacements_solved"] := true;

    # Add veils
    if keep_veiled then
      struct["veils"] := struct["veils"] union veils;
      veiling_idx     := veiling_idx + 1;
      veiling_label   := cat('_V', veiling_idx);
    end if;
  end if;

  PrintEndProc(procname);
  return struct;
end proc: # SolveStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

HyperstaticSolver := proc(
  objs::{ # Structural objects
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set( {BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  exts::{ # External actions
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  },
  vars::{list},       # Variables
  hyper_vars::{list}, # Hyperstatic variables
  hyper_disp::{list}, # Hyperstatic displacements
  {
    timoshenko_beam::{boolean} := false, # Timoshenko beam flag
    implicit::{boolean}        := false  # Implicit flag
  }, $)

  description "Solve hyperstatic structure with inputs objects <objs>, external "
    "actions <exts>, variables <vars>, hyperstatic variables <hyper_vars>, "
    "hyperstatic displacements <hyper_disp> and optional Timoshenko beam flag "
    "<timoshenko_beam>.";

  local hyper_eq, hyper_load, hyper_comps, hyper_compliant_disp, hyper_support,
    i, obj, iso_vars, iso_sol, iso_eq, hyper_sol, P_energy, S_objs, E_objs, sol;
  PrintStartProc(procname);

  # Parse input objects and find objects with internal actions property
  S_objs := [seq(
    `if`(IsBeam(objs[i]) or IsRod(objs[i]), objs[i], NULL),
    i = 1..nops(objs))
    ];

  # Create a solution as function of the hyperstatic variables
  iso_vars := [seq(
    `if`(member(vars[i], hyper_vars), NULL, vars[i]),
    i = 1..nops(vars))
    ];
  iso_sol, iso_eq, iso_vars := IsostaticSolver(objs, exts, iso_vars, parse("implicit") = implicit);

  # Compute internal actions
  ComputeInternalActions(S_objs, exts, iso_sol);

  # Extract the deformable objects
  E_objs := [seq(
    `if`(not IsRigidBody(objs[i]), objs[i], NULL),
    i = 1..nops(objs))
    ];

  # Compute structure internal energy
  P_energy := ComputePotentialEnergy(
    E_objs, iso_sol,
    parse("timoshenko_beam") = timoshenko_beam
  );
  P_energy := Simplify(P_energy);

  # Compute the hyperstatic equation
  if keep_veiled then
    hyper_eq := Diff~(P_energy, hyper_vars) =~ hyper_disp, parse("veils") = iso_sol[-1];
  else
    hyper_eq := diff~(P_energy, hyper_vars) =~ hyper_disp;
  end if;

  # Check for implicit solution flag
  if (implicit) then
    hyper_sol := iso_sol;
  else
    if (verbose_mode > 1) then
      printf("%*sMessage (in HyperstaticSolver) solving the hyperstatic variables...\n", print_indent, "|   ");
    end if;
    # Solve hyperstatic equations
    hyper_sol := op(RealDomain[solve](convert(hyper_eq, signum), hyper_vars));
    if (hyper_sol = NULL) then
      error "HyperstaticSolver: hyperstatic solution not found";
    end if;

    # Substitute hyper_sol in P_energy
    P_energy := subs(hyper_sol, P_energy);

    if (verbose_mode > 1) then
      printf("%*sDONE\n", print_indent, "|   ");
    end if;
    sol := hyper_sol union subs(hyper_sol, iso_sol);
  end if;

  PrintEndProc(procname);
  if (_nresults = 4) then
    return sol, iso_eq union hyper_eq, iso_vars union hyper_vars, P_energy;
  else
    return sol;
  end
end proc: # HyperstaticSolver

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputePotentialEnergy := proc(
  objs::{ # Structure objects
    list({BEAM, ROD, SUPPORT, JOINT}),
    set( {BEAM, ROD, SUPPORT, JOINT})
  },
  sol::{list, set} := [], # Substitutions
  {
    timoshenko_beam::{boolean} := false, # Timoshenko beam flag
    dummy_vars::{list, set}    := []     # Dummy variables
  }, $)

  description "Compute the internal potential energy of the structure given the "
    "objects <objs> and optional Timoshenko beam flag <timoshenko_beam>.";

  local dummy_vars_subs, obj, P, x, f, FJX, FJY, FJZ, MJX, MJY, MJZ;
  PrintStartProc(procname);

  dummy_vars_subs := dummy_vars =~ [seq(0, i = 1..nops(dummy_vars))];

  P := 0;
  for obj in objs do
    if IsBeam(obj) or IsRod(obj) then
      # Normal action N contribution
      if (member(N, map(lhs, obj["internal_actions"]))) and
          (subs(obj["internal_actions"](x), N(x)) <> 0) then
        subs(obj["internal_actions"](x), N(x));
        subs(dummy_vars_subs, %);
        P := P + integrate(
            eval(`if`(nops(dummy_vars) > 0, 2 * (%% - %) * %, %%^2))/
              (2*obj["material"]["elastic_modulus"]*obj["area"](x)),
            x = 0..obj["length"]);
      end if;
      if timoshenko_beam then
        # Shear action Ty contribution
        if (member(Ty, map(lhs, obj["internal_actions"]))) and
            (subs(obj["internal_actions"](x), Ty(x)) <> 0) then
          subs(obj["internal_actions"](x), Ty(x));
          subs(dummy_vars_subs, %);
          P := P + integrate(
              eval(`if`(nops(dummy_vars) > 0, 2 * (%% - %) * %, %%^2))/
                (2*obj["timo_shear_coeff"](x)[1]*obj["material"]["shear_modulus"]*obj["area"](x)),
              x = 0..obj["length"]);
        end if;
        # Shear action Tz contribution
        if (member(Tz, map(lhs, obj["internal_actions"]))) and
            (subs(obj["internal_actions"](x), Tz(x)) <> 0) then
          subs(obj["internal_actions"](x), Tz(x));
          subs(dummy_vars_subs, %);
          P := P + integrate(
              eval(`if`(nops(dummy_vars) > 0, 2 * (%% - %) * %, %%^2))/
                (2*obj["timo_shear_coeff"](x)[2]*obj["material"]["shear_modulus"]*obj["area"](x)),
              x = 0..obj["length"]);
        end if;
      end if;
      # Bending moment action Mx contribution
      if (member(Mx, map(lhs, obj["internal_actions"]))) and
          (subs(obj["internal_actions"](x), Mx(x)) <> 0) then
        subs(obj["internal_actions"](x), Mx(x));
        subs(dummy_vars_subs, %);
        P := P + integrate(
            eval(`if`(nops(dummy_vars) > 0, 2 * (%% - %) * %, %%^2))/
              (2*obj["material"]["shear_modulus"]*obj["inertias"][1](x)),
            x = 0..obj["length"]);
          end if;
      # Bending moment action My contribution
      if (member(My, map(lhs, obj["internal_actions"]))) and
          (subs(obj["internal_actions"](x), My(x)) <> 0) then
        subs(obj["internal_actions"](x), My(x));
        subs(dummy_vars_subs, %);
        P := P + integrate(
            eval(`if`(nops(dummy_vars) > 0, 2 * (%% - %) * %, %%^2))/
              (2*obj["material"]["elastic_modulus"]*obj["inertias"][2](x)),
            x = 0..obj["length"]);
      end if;
      # Bending moment action Mz contribution
      if (member(Mz, map(lhs, obj["internal_actions"]))) and
          (subs(obj["internal_actions"](x), Mz(x)) <> 0) then
        subs(obj["internal_actions"](x), Mz(x));
        subs(dummy_vars_subs, %);
        P := P + integrate(
            eval(`if`(nops(dummy_vars) > 0, 2 * (%% - %) * %, %%^2))/
              (2*obj["material"]["elastic_modulus"]*obj["inertias"][3](x)),
            x = 0..obj["length"]);
      end if;
    elif IsCompliantSupport(obj) then
      # Support reaction Fx contribution
      if (subs(obj["support_reactions"], FX) <> 0) and
          (obj["stiffness"](x)[1] <> infinity) and
          (obj["constrained_dof"][1] <> 0) then
        P := P + Subs(obj["support_reactions"], sol, ComputeSpringEnergy(-FX, (x -> obj["stiffness"](x)[1])));
      end if;
      # Support reaction Fy contribution
      if (subs(obj["support_reactions"], FY) <> 0) and
          (obj["stiffness"](x)[2] <> infinity) and
          (obj["constrained_dof"][2] <> 0) then
        P := P + Subs(obj["support_reactions"], sol, ComputeSpringEnergy(-FY, (x -> obj["stiffness"](x)[2])));
      end if;
      # Support reaction Fz contribution
      if (subs(obj["support_reactions"], FZ) <> 0) and
          (obj["stiffness"](x)[3] <> infinity) and
          (obj["constrained_dof"][3] <> 0) then
        P := P + Subs(obj["support_reactions"], sol, ComputeSpringEnergy(-FZ, (x -> obj["stiffness"](x)[3])));
      end if;
      # Support reaction Mx contribution
      if (subs(obj["support_reactions"], MX) <> 0) and
          (obj["stiffness"](x)[4] <> infinity) and
          (obj["constrained_dof"][4] <> 0) then
        P := P + Subs(obj["support_reactions"], sol, ComputeSpringEnergy(-MX, (x -> obj["stiffness"](x)[4])));
      end if;
      # Support reaction My contribution
      if (subs(obj["support_reactions"], MY) <> 0) and
          (obj["stiffness"](x)[5] <> infinity) and
          (obj["constrained_dof"][5] <> 0) then
        P := P + Subs(obj["support_reactions"], sol, ComputeSpringEnergy(-MY, (x -> obj["stiffness"](x)[5])));
      end if;
      # Support reaction Mz contribution
      if (subs(obj["support_reactions"], MZ) <> 0) and
          (obj["stiffness"](x)[6] <> infinity) and
          (obj["constrained_dof"][6] <> 0) then
        P := P + Subs(obj["support_reactions"], sol, ComputeSpringEnergy(-MZ, (x -> obj["stiffness"](x)[6])));
      end if;
    elif IsCompliantJoint(obj) then
      # Joint forces along X axis
      if (obj["stiffness"](x)[1] <> infinity) and
           (obj["constrained_dof"][1] <> 0) then
        FJX := 0;
        # Get all the forces along X axis
        for f in obj["forces"] do
          if member(f["target"], obj["shell_targets"]) then
            FJX := FJX + f["components"][1];
          end if;
        end do;
        P := P + Subs(sol, ComputeSpringEnergy(FJX, (x -> obj["stiffness"](x)[1])));
      end if;
      # Joint forces along Y axis
      if (obj["stiffness"](x)[2] <> infinity) and
           (obj["constrained_dof"][2] <> 0) then
        FJY := 0;
        # Get all the forces along Y axis
        for f in obj["forces"] do
          if member(f["target"], obj["shell_targets"]) then
            FJY := FJY + f["components"][2];
          end if;
        end do;
        P := P + Subs(sol, ComputeSpringEnergy(FJY, (x -> obj["stiffness"](x)[2])));
      end if;
      # Joint forces along Z axis
      if (obj["stiffness"](x)[3] <> infinity) and
           (obj["constrained_dof"][3] <> 0) then
        FJZ := 0;
        # Get all the forces along Z axis
        for f in obj["forces"] do
          if member(f["target"], obj["shell_targets"]) then
            FJZ := FJZ + f["components"][3];
          end if;
        end do;
        P := P + Subs(sol, ComputeSpringEnergy(FJZ, (x -> obj["stiffness"](x)[3])));
      end if;
      # Joint moments along X axis
      if (obj["stiffness"](x)[4] <> infinity) and
           (obj["constrained_dof"][4] <> 0) then
        MJX := 0;
        # Get all the moments along X axis
        for f in obj["moments"] do
          if member(f["target"], obj["shell_targets"]) then
            MJX := MJX + f["components"][1];
          end if;
        end do;
        P := P + Subs(sol, ComputeSpringEnergy(MJX, (x -> obj["stiffness"](x)[4])));
      end if;
      # Joint moments along Y axis
      if (obj["stiffness"](x)[5] <> infinity) and
           (obj["constrained_dof"][5] <> 0) then
        MJY := 0;
        # Get all the moments along Y axis
        for f in obj["moments"] do
          if member(f["target"], obj["shell_targets"]) then
            MJY := MJY + f["components"][2];
          end if;
        end do;
        P := P + Subs(sol, ComputeSpringEnergy(MJY, (x -> obj["stiffness"](x)[5])));
      end if;
      # Joint moments along Z axis
      if (obj["stiffness"](x)[6] <> infinity) and
           (obj["constrained_dof"][6] <> 0) then
        MJZ := 0;
        # Get all the moments along Z axis
        for f in obj["moments"] do
          if member(f["target"], obj["shell_targets"]) then
            MJZ := MJZ + f["components"][3];
          end if;
        end do;
        P := P + Subs(sol, ComputeSpringEnergy(MJZ, (x -> obj["stiffness"](x)[6])));
      end if;
    end if;
  end do;

  PrintEndProc(procname);
  return P;
end proc: # PotentialEnergy

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsostaticSolver := proc(
  objs::{ # Structural objects
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set( {BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  exts::{ # External actions
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  },
  vars::{list}, # Variables to solve
  {
    implicit::{boolean} := false  # Implicit solution
  },
  $)

  description "Solve the isostatic structure equilibrium equation system given "
    "the structure objects <objs>, the external actions <exts> and the variables "
    "<vars> to solve.";

  local iso_eq, iso_eq_tmp, exts_comps, x, i, j, active_ext, iso_sol, A, B, rank_eq, iso_vars, PivotStrategy;
  PrintStartProc(procname);

  # Compute structure equations
  if (verbose_mode > 1) then
    printf("%*sMessage (in IsostaticSolver) computing the equilibrium equation for "
      "the isostatic structure...\n", print_indent, "|   ");
  end if;
  iso_eq := [];
  for i from 1 to nops(objs) do
    active_ext := {};
    for j from 1 to nops(exts) do
      if (exts[j]["target"] = objs[i]["name"]) then
        active_ext := active_ext union {exts[j]};
      end if;
    end do;
    iso_eq := iso_eq union NewtonEuler(active_ext, objs[i]);
    # Add joints and supports constraint equations
    if IsSupport(objs[i]) or IsJoint(objs[i]) then
      iso_eq := iso_eq union objs[i]["constraint_loads"];
    end if;
  end do;

  # Remove NULL equations
  iso_eq := remove(x -> x = 0, Simplify(iso_eq));

  # Remove equation related to rigid body motions
  # FIXME: for qloads should be different (this check can be skipped)
  iso_eq_tmp := iso_eq;
  exts_comps := map(x -> op(x["components"]), GetObjsByType({FORCE, MOMENT}, exts));
  #iso_eq := remove(x -> (member(0., (abs~(exts_comps) -~ abs(x)) *~ 1.) and (not member(0., (abs~(vars) -~ abs(x)) *~ 1.))), iso_eq_tmp);
  iso_eq := remove(x -> (member(0., (abs~(exts_comps) -~ abs(x)) *~ 1.) and (not has(vars, indets(x)))), iso_eq_tmp);
  if (not suppress_warnings) then
    if (nops(iso_eq_tmp) <> nops(iso_eq)) then
      WARNING("Message (in IsostaticSolver) the following list of equations were "
        "removed because they are related to rigid body motions:\n%1",
        convert(iso_eq_tmp *~ 1., set) minus convert(iso_eq *~ 1., set));
    end if;
  end if;

  # Remove non used variables
  iso_vars := remove(x -> (not has(iso_eq, x)), vars);
  if (not suppress_warnings) then
    if nops(iso_vars) <> nops(vars) then
      WARNING(
        "Message (in IsostaticSolver) the following list of variables were removed "
        "because they are not used in the structure equilibrium equations:\n%1",
        convert(vars, set) minus convert(iso_vars, set));
    end
  end if;

  if (verbose_mode > 1) then
    printf("%*sDONE\n", print_indent, "|   ");
  end if;

  if (verbose_mode > 1) then
    printf("%*sMessage (in IsostaticSolver) structure equilibrium equations:\n", print_indent, "|   ");
    print(<op(iso_eq)>);
    printf("%*sMessage (in IsostaticSolver) structure unknown variables:\n", print_indent, "|   ");
    print(iso_vars);
  end if;

  # Check for implicit solution flag
  if (implicit) then
    iso_sol := [];
  else
    # Matrix form
    A, B := LinearAlgebra[GenerateMatrix](iso_eq, iso_vars);
    A := Matrix(A, storage = sparse);

    if nops(iso_eq) = nops(iso_vars) then

      if (verbose_mode > 1) then
        printf("%*sMessage (in IsostaticSolver) A matrix visualization of the linear system:\n", print_indent, "|   ");
        print(plots[sparsematrixplot](A, matrixview));
      end if;
      if (verbose_mode > 1) then
        printf("%*sMessage (in IsostaticSolver) computing the structure reaction forces...\n", print_indent, "|   ");
      end if;
      # Solve structure equations (LinearSolver)
      #iso_sol := LinearSolver(iso_eq, iso_vars);
      iso_sol := op(RealDomain[solve](iso_eq, iso_vars)); # FIXME: this is a temporary fix (use LULEM when stable)
    else
      if (not suppress_warnings) then
        WARNING("Message (in IsostaticSolver) the system of equations is not "
          "consistent, trying to solve the system of equations anyway without "
          "LinearSolver");
      end if;
      # Solve structure equations (solve)
      iso_sol := op(RealDomain[solve](iso_eq, iso_vars));
      # append an empty list to the solution if the keep_veiled flag is set to true
      if keep_veiled then
        iso_sol := iso_sol union [[]];
      end if;
    end if;

    if iso_sol = NULL then
      error "isostaticSolver: isostatic solution not found";
    end if;

    if (verbose_mode > 1) then
      printf("%*sDONE\n", print_indent, "|   ");
    end if;
  end if;

  PrintEndProc(procname);
  if _nresults = 3 then
    return iso_sol, iso_eq, iso_vars;
  else
    return iso_sol;
  end
end proc: # IsostaticSolver

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeInternalActions := proc(
  objs::{ # Structure objects
    list({BEAM, ROD}),
    set( {BEAM, ROD})
  },
  exts::{ # External actions
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  },
  sol::{list, set}, # Structure solution
  $)::{nothing};

  description "Programmatic computation of internal actions for structure"
    "objects with given external actions and structure solution.";

  local i, j, active_ext, subs_ext;
  PrintStartProc(procname);

  # Substitute structure solution into loads
  subs_ext := map(convert, map2(Subs, sol, map(op, exts)), table);

  for i from 1 to nops(objs) do
    # Extract active loads
    active_ext := {};
    for j from 1 to nops(subs_ext) do
      if (subs_ext[j]["target"] = objs[i]["name"]) then
        active_ext := active_ext union {subs_ext[j]};
      end if;
    end do;
    # Compute internal actions
    InternalActions(objs[i], active_ext);
  end do;

  PrintEndProc(procname);
end proc: # ComputeInternalActions

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

InternalActions := proc(
  obj::{BEAM, ROD}, # Structure object
  exts::{           # External actions
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  },
  $)::{nothing};

  description "Programmatic computation of internal actions for structure "
    "object with given external actions, it returns the internal actions as "
    "function of the axial variable 'x'.";

  local i, ia, N_sol, Ty_sol, Tz_sol, Mx_sol, My_sol, Mz_sol, x;
  PrintStartProc(procname);

  N_sol  := 0;
  Ty_sol := 0;
  Tz_sol := 0;
  Mx_sol := 0;
  My_sol := 0;
  Mz_sol := 0;

  # Assumptions
  # NOTE: assumptions higly help readability of the solution and improve
  # computation time, but results must be considered valid only in the assumed range
  Physics[Assume](x > 0, x < obj["length"]);

  # Compute internal actions for concentrated loads as effect overlay
  for i from 1 to nops(exts) do
    if IsForce(exts[i]) then
      N_sol  := Simplify(N_sol  - piecewise(x >= exts[i]["coordinate"][1] and x <= obj["length"], exts[i]["components"][1]), piecewise);
      Ty_sol := Simplify(Ty_sol + piecewise(x >= exts[i]["coordinate"][1] and x <= obj["length"], exts[i]["components"][2]), piecewise);
      Tz_sol := Simplify(Tz_sol + piecewise(x >= exts[i]["coordinate"][1] and x <= obj["length"], exts[i]["components"][3]), piecewise);
      My_sol := Simplify(My_sol - integrate(piecewise(x >= exts[i]["coordinate"][1] and x <= obj["length"], exts[i]["components"][3]), x = 0..x), piecewise);
      Mz_sol := Simplify(Mz_sol + integrate(piecewise(x >= exts[i]["coordinate"][1] and x <= obj["length"], exts[i]["components"][2]), x = 0..x), piecewise);
    elif IsMoment(exts[i]) then
      Mx_sol := Simplify(Mx_sol - piecewise(x >= exts[i]["coordinate"][1] and x <= obj["length"], exts[i]["components"][1]), piecewise);
      My_sol := Simplify(My_sol - piecewise(x >= exts[i]["coordinate"][1] and x <= obj["length"], exts[i]["components"][2]), piecewise);
      Mz_sol := Simplify(Mz_sol - piecewise(x >= exts[i]["coordinate"][1] and x <= obj["length"], exts[i]["components"][3]), piecewise);
    elif IsQForce(exts[i]) then
      N_sol  := Simplify(N_sol  - integrate(exts[i]["components"](x)[1], x = 0..x), piecewise);
      Ty_sol := Simplify(Ty_sol + integrate(exts[i]["components"](x)[2], x = 0..x), piecewise);
      Tz_sol := Simplify(Tz_sol + integrate(exts[i]["components"](x)[3], x = 0..x), piecewise);
      My_sol := Simplify(My_sol - integrate(integrate(exts[i]["components"](x)[3], x = 0..x), x = 0..x), piecewise);
      Mz_sol := Simplify(Mz_sol + integrate(integrate(exts[i]["components"](x)[2], x = 0..x), x = 0..x), piecewise);
    elif IsQMoment(FMQ[i]) then
      Mx_sol := Simplify(Mx_sol - integrate(exts[i]["components"](x)[1], x = 0..x), piecewise);
      My_sol := Simplify(My_sol - integrate(exts[i]["components"](x)[2], x = 0..x), piecewise);
      Mz_sol := Simplify(Mz_sol - integrate(exts[i]["components"](x)[3], x = 0..x), piecewise);
    end if;
  end do;

  ia := [
    N  = unapply( N_sol, x), Ty = unapply(Ty_sol, x), Tz = unapply(Tz_sol, x),
    Mx = unapply(Mx_sol, x), My = unapply(My_sol, x), Mz = unapply(Mz_sol, x)
    ];

  if IsRod(obj) then
    ia := [ia[1]];
  end if;

  if (verbose_mode > 1) then
  printf(
    "%*sMessage (in InternalActions) updating %s %s's internal actions...\n",
    print_indent, "|   ", obj["type"], obj["name"]
    );
  end if;

  obj["internal_actions"] := ia;

  if (verbose_mode > 1) then
    printf("%*sDONE\n", print_indent, "|   ");
  end if;

  PrintEndProc(procname);
end proc: # InternalActions

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeDisplacements := proc(
  objs::{ # Structural objects
    list({BEAM, ROD, SUPPORT, JOINT}),
    set( {BEAM, ROD, SUPPORT, JOINT})
  },
  exts::{ # External loads
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  },
  sol::{list, set}, # Solution of the structure
  {
    timoshenko_beam::{boolean} := false # Timoshenko beam flag
  }, $)::{nothing};

  description "Compute the structure displacements and rotations.";

  local obj, x, X, disp, rx_sol, ry_sol, rz_sol, ux_sol, uy_sol, uz_sol;
  PrintStartProc(procname);

  # Cicle on the structure objects
  for obj in objs do
    # Beam
    if IsBeam(obj) then
      #Physics[Assume](X > 0, X < obj["length"]);
      # Compute displacements
      rx_sol :=  integrate(subs(obj["internal_actions"](x), Mx(x)/(obj["material"]["shear_modulus"]*obj["inertias"][1](x))), x = 0..X);
      ry_sol :=  integrate(subs(obj["internal_actions"](x), My(x)/(obj["material"]["elastic_modulus"]*obj["inertias"][2](x))), x = 0..X);
      rz_sol :=  integrate(subs(obj["internal_actions"](x), Mz(x)/(obj["material"]["elastic_modulus"]*obj["inertias"][3](x))), x = 0..X);
      ux_sol :=  integrate(subs(obj["internal_actions"](x), N(x)/(obj["material"]["elastic_modulus"]*obj["area"](x))), x = 0..X);
      uy_sol :=  integrate(eval(rz_sol, X = x), x = 0..X);
      uz_sol := -integrate(eval(ry_sol, X = x), x = 0..X);
      if timoshenko_beam then
        uy_sol := uy_sol + integrate(subs(obj["internal_actions"](x), Ty(x)/(obj["timo_shear_coeff"](x)[1]*obj["material"]["shear_modulus"]*obj["area"](x))), x = 0..X);
        uz_sol := uz_sol + integrate(subs(obj["internal_actions"](x), Tz(x)/(obj["timo_shear_coeff"](x)[2]*obj["material"]["shear_modulus"]*obj["area"](x))), x = 0..X);
      end if;
      disp := [
        ux = unapply(ux_sol, X), uy = unapply(uy_sol, X), uz = unapply(uz_sol, X),
        rx = unapply(rx_sol, X), ry = unapply(ry_sol, X), rz = unapply(rz_sol, X)
      ];

      # Update object displacements
      obj["displacements"] := disp;

    # Rod
    elif IsRod(obj) then
      #Physics[Assume](X > 0, X < obj["length"]);
      # Compute displacements
      ux_sol := integrate(subs(obj["internal_actions"](x), N(x)/(obj["material"]["elastic_modulus"]*obj["area"](x))), x = 0..X);
      disp := [
        ux = unapply(ux_sol, X)
      ];

      # Update object displacements
      obj["displacements"] := disp;

    # Support
    elif IsCompliantSupport(obj) then
      # Compute displacements
      ComputeSupportDisplacements(obj);

    # Joint
    elif IsCompliantJoint(obj) then
      # Compute displacements
      ComputeJointDisplacements(obj, sol);
    end if;
  end do;

  PrintEndProc(procname);
end proc: # ComputeDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputePunctualDisplacement := proc(
  struct::{STRUCTURE}, # Structure
  objs::{ # Object on which the coordinates are defined
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  coords::{list},     # Punctual coordinates defined in obj reference frame
  directions::{list}, # Displacement directions defined in obj reference frame
  RFs::{list},        # Reference frames for the directions
  {
    timoshenko_beam::{boolean} := false, # Timoshenko beam flag
    unveil_results::{boolean}  := true   # Unveil results flag
  },
  $)::{list};

  description "Compute the Structure <struct> punctual displacements of the "
    "object <obj> at the coordinates <coords> in the directions <directions>. "
    "The directions are defined in the reference frame <RFs>. Optional argument "
    "<timoshenko_beam> is a boolean flag to use Timoshenko beam model, <unveil_results> "
    "is a boolean flag to unveil the results.";

  local out, struct_copy, obj, objs_names, dummy_loads, subs_obj, obj_coords,
    obj_targets, x, subs_null_dummy, disp, i, j, d_coords, sw_tmp;
  PrintStartProc(procname);

  # Substitute -1 entries of coords with the corresponding object length
  d_coords := [seq(`if`(coords[i] = -1, objs[i]["length"], coords[i]), i = 1..nops(coords))];

  # Set module local variable keep_veiled
  keep_veiled := not unveil_results;

  # Create a copy of the structure
  struct_copy := CopyStructure(struct);

  # Get objects names
  objs_names := GetNames(objs);

  # Disable warnings temporarily FIXME: this is a temporary fix
  sw_tmp := suppress_warnings;
  suppress_warnings := true;

  # Replace Rods with Beams to be able to compute the displacements in all directions
  for obj in map(GetObjByName, objs_names, struct_copy["objects"]) do
    if IsRod(obj) then
      subs_obj := MakeBeam(obj["name"], obj["length"], obj["frame"], parse("area") = obj["area"], parse("material") = obj["material"]);
      # Remove load on unconstrained direction
      subs_obj["admissible_loads"] := [1, 1, 1, 0, 1, 1];
      # Replace object in struct_copy
      struct_copy["objects"] := remove(x -> x["name"] = obj["name"], struct_copy["objects"]);
      struct_copy["objects"] := struct_copy["objects"] union {eval(subs_obj)};
    end if;
  end do;

  # Update struct_copy supports and joints for the new objects
  for obj in struct_copy["objects"] do
    if IsSupport(obj) then
      # Re-Make support to generate new loads and constraint compliant with substituted objects (joint is made because earth is already in the list of targets)
      obj_targets := map(GetObjByName, remove(x -> x = earth["name"], obj["targets"]), struct_copy["objects"]);
      obj_coords  := obj["coordinates"][2..-1];
      subs_obj := MakeSupport(obj["name"], obj["constrained_dof"], obj_targets, obj_coords, obj["frame"], parse("stiffness") = obj["stiffness"]);
      # Replace object in struct_copy
      struct_copy["objects"] := remove(x -> x["name"] = obj["name"], struct_copy["objects"]);
      struct_copy["objects"] := struct_copy["objects"] union {eval(subs_obj)};
    elif IsJoint(obj) then
      # Re-Make joint to generate new loads and constraint compliant with substituted objects
      obj_targets := map(GetObjByName, obj["targets"], struct_copy["objects"]);
      subs_obj := MakeJoint(obj["name"], obj["constrained_dof"], obj_targets, obj["coordinates"], obj["frame"], parse("stiffness") = obj["stiffness"]);
      # Replace object in struct_copy
      struct_copy["objects"] := remove(x -> x["name"] = obj["name"], struct_copy["objects"]);
      struct_copy["objects"] := struct_copy["objects"] union {eval(subs_obj)};
    end if;
  end do;

  # Create dummy loads in the directions of interest
  subs_null_dummy := [];
  for i from 1 to nops(objs_names) do
    dummy_loads := eval~([
    `if`(directions[i,1] = 1, MakeForce( [dFx_||i,0,0], d_coords[i], GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL),
    `if`(directions[i,2] = 1, MakeForce( [0,dFy_||i,0], d_coords[i], GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL),
    `if`(directions[i,3] = 1, MakeForce( [0,0,dFz_||i], d_coords[i], GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL),
    `if`(directions[i,4] = 1, MakeMoment([dMx_||i,0,0], d_coords[i], GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL),
    `if`(directions[i,5] = 1, MakeMoment([0,dMy_||i,0], d_coords[i], GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL),
    `if`(directions[i,6] = 1, MakeMoment([0,0,dMz_||i], d_coords[i], GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL)
    ]);

    # Null dummy loads substitution list
    subs_null_dummy := subs_null_dummy union ([dFx_||i, dFy_||i, dFz_||i, dMx_||i, dMy_||i, dMz_||i] =~ [0, 0, 0, 0, 0, 0]);

    # Add dummy loads to the structure copy
    struct_copy["external_actions"] := struct_copy["external_actions"] union dummy_loads;
  end do;

  StoredData := StoredData union subs_null_dummy;

  # Solve the structure copy
  SolveStructure(
    struct_copy,
    parse("compute_internal_actions") = false,
    parse("compute_displacements")    = false,
    parse("compute_potential_energy") = true,
    parse("timoshenko_beam")          = timoshenko_beam,
    parse("implicit")                 = false,
    parse("unveil_results")           = unveil_results,
    parse("dummy_vars")               = lhs~(subs_null_dummy)
    );

  # Enable warnings
  suppress_warnings := sw_tmp;

  # Compute punctual displacements
  out := [];
  for i from 1 to nops(objs_names) do
    disp := [];
    if directions[i,1] = 1 then
      disp := disp union [ux = subs(subs_null_dummy, Diff(struct_copy["potential_energy"], dFx_||i, parse("veils") = struct_copy["veils"]))];
    end if;
    if directions[i,2] = 1 then
      disp := disp union [uy = subs(subs_null_dummy, Diff(struct_copy["potential_energy"], dFy_||i, parse("veils") = struct_copy["veils"]))];
    end if;
    if directions[i,3] = 1 then
      disp := disp union [uz = subs(subs_null_dummy, Diff(struct_copy["potential_energy"], dFz_||i, parse("veils") = struct_copy["veils"]))];
    end if;
    if directions[i,4] = 1 then
      disp := disp union [rx = subs(subs_null_dummy, Diff(struct_copy["potential_energy"], dMx_||i, parse("veils") = struct_copy["veils"]))];
    end if;
    if directions[i,5] = 1 then
      disp := disp union [ry = subs(subs_null_dummy, Diff(struct_copy["potential_energy"], dMy_||i, parse("veils") = struct_copy["veils"]))];
    end if;
    if directions[i,6] = 1 then
      disp := disp union [rz = subs(subs_null_dummy, Diff(struct_copy["potential_energy"], dMz_||i, parse("veils") = struct_copy["veils"]))];
    end if;
    # Select displacement relative to desired directions and add to the list of displacements
    out := out union [disp];
  end do;

   struct_copy["veils"] := subs(subs_null_dummy,  struct_copy["veils"]);

  # Simplify output
  out := Simplify(out);

  PrintEndProc(procname);
  if _nresults = 1 then
    return out;
  else
    return out, struct_copy["veils"];
  end if;
end proc: # ComputePunctualDisplacement

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeObjectFrameDisplacements := proc(
  struct::STRUCTURE, # Structure to compute the total displacements
  {
    timoshenko_beam::{boolean} := false, # Timoshenko beam flag
    unveil_results::{boolean}  := true   # Unveil results flag
  },
  $)::{STRUCTRURE};

description "Compute the total displacements of the structure <struct>.";

  local i, RF_nt, nx, ny, nz, theta, subs_n, subs_t, disp, veils, objs;
  PrintStartProc(procname);

  # Compute punctual displacements at origin for all the structure objects
  # FIXME: this is overkill, we should compute only the displacements in the
  # direction of joints dof (and not even all of them)
  disp, veils := ComputePunctualDisplacement(struct,
                                             convert(struct["objects"],list),
                                             [seq([0, 0, 0], i = 1..nops(struct["objects"]))],
                                             [seq([1, 1, 1, 1, 1, 1], i = 1..nops(struct["objects"]))],
                                             map(x-> x["frame"], convert(struct["objects"], list)),
                                             parse("timoshenko_beam") = timoshenko_beam,
                                             parse("unveil_results")  = unveil_results
                                             ):

  # Update objects <frame_displacements>
  for i from 1 to nops(struct["objects"]) do
    # FIXME: not all displacement are necessarily computed
    # Rotation about generic axis
    RF_nt  := Matrix(4, 4, [[-nx^2*cos(theta) + nx^2 + cos(theta), -nx*ny*cos(theta) - sin(theta)*nz + nx*ny, -nx*nz*cos(theta) + sin(theta)*ny + nx*nz, 0], [-nx*ny*cos(theta) + sin(theta)*nz + nx*ny, (-nx^2 + nz^2 + 1)*cos(theta) + nx^2 - nz^2, -sin(theta)*nx - ny*nz*(cos(theta) - 1), 0], [-nx*nz*cos(theta) - sin(theta)*ny + nx*nz, sin(theta)*nx - ny*nz*(cos(theta) - 1), -cos(theta)*(-2*nx^2 + nz^2) - 2*nx^2 + nz^2 + 1, 0], [0, 0, 0, 1]]);
    if subs(disp[i],Norm2([rx, ry, rz])) <> 0. and subs(disp[i],Norm2([rx, ry, rz])) <> 0 then
      subs_n := [nx, ny, nz] =~ subs(disp[i], [rx, ry, rz] /~ Norm2([rx, ry, rz]));
    else
      subs_n := [nx, ny, nz] =~ [0, 0, 1];
    end if;
    subs_t := theta = subs(disp[i], Norm2([rx, ry, rz]));
    struct["objects"][i]["frame_displacements"] := Translate(op(subs(disp[i], [ux, uy, uz]))).
                                                   subs(subs_n, subs_t, RF_nt);
  end do;

  PrintEndProc(procname);
  return veils;
end proc; # ComputeObjectFrameDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CopyStructure := proc(
  struct::{STRUCTURE}, # Structure to copy
  $)::{STRUCTURE};

description "Create a copy of the structure <struct> and its objects.";

  local struct_copy, obj, action;
  PrintStartProc(procname);

  # Create a copy of the structure
  struct_copy := copy(struct);

  # Substitute objects in the structure with a copy
  for obj in struct_copy["objects"] do
    struct_copy["objects"] := remove(x -> x["name"] = obj["name"], struct_copy["objects"]);
    struct_copy["objects"] := struct_copy["objects"] union {copy(obj)};
  end do;

  # Substitute external actions in the structure with a copy
  for action in struct_copy["external_actions"] do
    struct_copy["external_actions"] := remove(x -> x["name"] = action["name"], struct_copy["external_actions"]);
    struct_copy["external_actions"] := struct_copy["external_actions"] union {copy(action)};
  end do;

  PrintEndProc(procname);
  return struct_copy;
end proc: # CopyStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

LinearSolver := proc(
  eqns::{list, set}, # Equations
  vars::{list, set}, # Variables
  $)::{list};

  description "Solve the linear system of equations <eqns> for the variables "
  "<vars>.";

  local T, sol, sol_tmp, A, b, _Q, PivotingStrategy;
  PrintStartProc(procname);

  # Matrix form of the linear system
  A, b := LinearAlgebra[GenerateMatrix](eqns, vars);

  use LULEM in
  LULEM[SetVerbosity](true);
  #LULEM[AssignData](StoredData);
  if has(map(type, A, constant), false) then
    PivotingStrategy := PivotingStrategy_Slength;
  else
    PivotingStrategy := PivotingStrategy_numeric;
  end if;
  T := LU(A, veiling_label);
  printf("Solved with rank: %d\n", T["rank"]);
  sol_tmp := SolveLinearSystem(T, b, veiling_label);
  if keep_veiled then
    # Remove indexed type from veils
    LEM[VeilList](veiling_label);
    lhs~(%) =~ map2(op,0,lhs~(%)) ||~ __ ||~ (op~(lhs~(%)));
    # Add veils to solution
    sol := convert(vars =~ subs(%, sol_tmp), list) union [subs(%,%%)];
  else
    sol := convert(vars =~ LEM[VeilSubs](sol_tmp, veiling_label), list);
  end if;
  ForgetVeil(veiling_label);
  #LULEM[ForgetData]();
  end use;

  PrintEndProc(procname);
  return Simplify(sol);
end proc: # LinearSolver

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ObjectColor := proc(
  obj::{BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH}, # Object to be colored
  $)::{string};

description "Return the color of the object <obj>.";

  local color;
  PrintStartProc(procname);

  if IsBeam(obj) then
    color := Beam_color;
  elif IsRod(obj) then
    color := Rod_color;
  elif IsRigidBody(obj) then
    color := RigidBody_color;
  elif IsCompliantSupport(obj) then
    color := CompliantSupport_color;
  elif IsSupport(obj) then
    color := Support_color;
  elif IsCompliantJoint(obj) then
    color := CompliantJoint_color;
  elif IsJoint(obj) then
    color := Joint_color;
  elif IsEarth(obj) then
    color := Earth_color;
  end if;

  PrintEndProc(procname);
  return color;
end proc: # ObjectColor


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotRigidBody := proc(
  obj::RIGID_BODY, # Rigid body to be plot
  joints::{ # Joint and support objects
    list({SUPPORT, JOINT}),
    set({SUPPORT, JOINT})
  },
  c_loads::{ # Concentrated loads
    list({FORCE, MOMENT}),
    set({FORCE, MOMENT})
  },
  {
    data::{list(`=`),set(`=`)} := [] # Substitutions
  },
  $)::{procedure};

  description "Plot the RIGID_BODY object <obj>.";

  local P1, P2, js, idx, lines, load, out;
  PrintStartProc(procname);

  lines := [];
  P1 := subs(op(data), Project([op(obj["COM"]), 1], obj["frame"], ground));
  for js in joints do
    member(obj["name"], js["targets"], 'idx');
    P2 := subs(op(data), Project([op(js["coordinates"][idx]), 1], obj["frame"], ground));
    lines := lines union [plottools:-line(convert(P1[1..3], list), convert(P2[1..3], list), thickness = 6)];
  end do;

  for load in c_loads do
    P2 :=  subs(op(data), Project([op(load["coordinate"]), 1], obj["frame"], ground));
    lines := lines union [plottools:-line(convert(P1[1..3], list), convert(P2[1..3], list), thickness = 6)];
  end do;

  out := plots:-display(lines, linestyle = solid, color = ObjectColor(obj), scaling = constrained);

  PrintEndProc(procname);
  return out;
end proc: # PlotRigidBody

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotDeformedRigidBody := proc(
  obj::RIGID_BODY, # Rigid body to be plot
  joints::{ # Joint and support objects
    list({SUPPORT, JOINT}),
    set({SUPPORT, JOINT})
  },
  c_loads::{ # Concentrated loads
    list({FORCE, MOMENT}),
    set({FORCE, MOMENT})
  },
  {
    data::{list(`=`),set(`=`)} := [], # Substitutions
    scaling::{numeric}         := 1.0 # Scaling factor
  },
  $)::{procedure};

  description "Plot the deformed RIGID_BODY object <obj>.";

  local P1, P2, rfd, js, idx, lines, load, out;
  PrintStartProc(procname);

  lines := [];
  rfd := subs(op(data), obj["frame"] . ((obj["frame_displacements"] - LinearAlgebra[IdentityMatrix](4)) *~ scaling + LinearAlgebra[IdentityMatrix](4)));
  P1 := subs(op(data), Project([op(obj["COM"]), 1], rfd, ground));
  for js in joints do
    member(obj["name"], js["targets"], 'idx');
    P2 := subs(op(data), Project([op(js["coordinates"][idx]), 1], rfd, ground));
    lines := lines union [plottools:-line(convert(P1[1..3], list), convert(P2[1..3], list), thickness = 6)];
  end do;

  for load in c_loads do
    P2 :=  subs(op(data), Project([op(load["coordinate"]), 1], rfd, ground));
    lines := lines union [plottools:-line(convert(P1[1..3], list), convert(P2[1..3], list), thickness = 6)];
  end do;

  out := plots:-display(lines, linestyle = solid, color = ObjectColor(obj), parse("scaling") = constrained);

  PrintEndProc(procname);
  return out;
end proc: # PlotDeformedRigidBody

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotBeam := proc(
  obj::{BEAM}, # Beam to be plot
  {
    data::{list(`=`),set(`=`)} := [] # Substitutions
  },
  $)::{procedure};

  description "Plot the BEAM object <obj>.";

  local P1, P2, out;
  PrintStartProc(procname);

  P1 := subs(op(data), Origin(obj["frame"]));
  P2 := subs(op(data), Project([obj["length"], 0, 0, 1], obj["frame"], ground));

  out := plots:-display(
    plottools:-line(convert(P1[1..3], list), convert(P2[1..3], list), thickness = 6),
    linestyle = solid, color = ObjectColor(obj), scaling = constrained);

  PrintEndProc(procname);
  return out;
end proc: # PlotBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotDeformedBeam := proc(
  obj::{BEAM}, # Beam to be plot
  {
    data::{list(`=`),set(`=`)} := [], # Substitutions
    scaling::{numeric}         := 1   # Scaling factor
  },
  $)::{procedure};

  description "Plot the BEAM object <obj>.";

  local sc, rfd, x, out;
  PrintStartProc(procname);

  rfd := subs(op(data), obj["frame"] . ((obj["frame_displacements"] - LinearAlgebra[IdentityMatrix](4)) *~ scaling + LinearAlgebra[IdentityMatrix](4)));

  sc := subs(op(data), Project(subs(obj["displacements"](x), [x, 0, 0, 0] +~ [ux(x) *~ scaling, uy(x) *~ scaling, uz(x) *~ scaling, 1]), rfd, ground)[1..3]);

  out := plots:-display(
    plots:-spacecurve(sc, x = subs(op(data), 0..obj["length"]), thickness = 6),
    linestyle = solid, color = ObjectColor(obj), parse("scaling") = constrained);

  PrintEndProc(procname);
  return out;
end proc: # PlotDeformedBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotRod := proc(
  obj::{ROD}, # Rod to be plot
  {
    data::{list(`=`),set(`=`)} := [] # Substitutions
  },
  $)::{procedure};

  description "Plot the ROD object <obj>.";

  local P1, P2, out;
  PrintStartProc(procname);

  P1 := subs(op(data), Origin(obj["frame"]));
  P2 := subs(op(data), Project([obj["length"], 0, 0, 1], obj["frame"], ground));

  out := plots:-display(
    plottools:-line(convert(P1[1..3], list), convert(P2[1..3], list), thickness = 4),
    linestyle = solid, color = ObjectColor(obj), scaling = constrained);

  PrintEndProc(procname);
  return out;
end proc: # PlotRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotDeformedRod := proc(
  obj::{ROD}, # Rod to be plot
  {
    data::{list(`=`),set(`=`)} := [], # Substitutions
    scaling::{numeric}         := 1   # Scaling factor
  },
  $)::{procedure};

  description "Plot the ROD object <obj>.";

  local P1, P2, rfd, out;
  PrintStartProc(procname);

  rfd := obj["frame"] . ((obj["frame_displacements"] - LinearAlgebra[IdentityMatrix](4)) *~ scaling + LinearAlgebra[IdentityMatrix](4));

  P1 := subs(op(data), Origin(rfd));
  P2 := subs(op(data), Project([obj["length"] + subs(obj["displacements"](obj["length"]), ux(obj["length"]) *~ scaling), 0, 0, 1], rfd, ground));

  out := plots:-display(
    plottools:-line(convert(P1[1..3], list), convert(P2[1..3], list), thickness = 4),
    linestyle = solid, color = ObjectColor(obj), parse("scaling") = constrained);

  PrintEndProc(procname);
  return out;
end proc: # PlotDeformedRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotJoint := proc(
  obj::{JOINT}, # Joint to be plot
  targets::{ # Joint targets
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  {
    data::{list(`=`),set(`=`)} := [] # Substitutions
  },
  $)::{procedure};

  description "Plot the JOINT object <obj>.";

  local O, out;
  PrintStartProc(procname);

  O := subs(op(data), Origin(
    GetObjByName(obj["targets"][1], targets)["frame"].
    Translate(op(ListPadding(obj["coordinates"][1],3)))
    ));

  out := plots:-display(
    plottools:-point(convert(O[1..3], list), symbol='solidsphere', symbolsize = 20),
    linestyle = solid, color = ObjectColor(obj), scaling = constrained);

  PrintEndProc(procname);
  return out;
end proc: # PlotJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotDeformedJoint := proc(
  obj::{JOINT}, # Joint to be plot
  targets::{ # Joint targets
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  {
    data::{list(`=`),set(`=`)} := [], # Substitutions
    scaling::{numeric}         := 1   # Scaling factor
  },
  $)::{procedure};

  description "Plot the JOINT object <obj>.";

  local O, rfd, out;
  PrintStartProc(procname);
  #TODO: add compliant joint deformation

  rfd := ((obj["frame_displacements"] - LinearAlgebra[IdentityMatrix](4)) *~ scaling + LinearAlgebra[IdentityMatrix](4));

  # FIXME: joint target may be a joint itself
  O := subs(op(data), Origin(
     GetObjByName(obj["targets"][1], targets)["frame"].
     Translate(op(ListPadding(obj["coordinates"][1],3))))[1..3] +~
     Project(Origin(rfd)[1..3], obj["frame"], ground)
    );

  out := plots:-display(
    plottools:-point(convert(O, list), symbol='solidsphere', symbolsize = 20),
    linestyle = solid, color = ObjectColor(obj), parse("scaling") = constrained);

  PrintEndProc(procname);
  return out;
end proc: # PlotDeformedJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotSupport := proc(
  obj::{SUPPORT}, # Support to be plotted
  targets::{ # Support targets
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  {
    data::{list(`=`),set(`=`)} := [] # Substitutions
  },
  $)::{procedure};

  description "Plot the SUPPORT object <obj>.";

  local O, out;
  PrintStartProc(procname);

  if (nops(obj["targets"]) > 1) then
    O := subs(op(data), Origin(
      GetObjByName(obj["targets"][2], targets)["frame"].
      Translate(op(ListPadding(obj["coordinates"][2],3)))
      ));
  else
    O := subs(op(data), Origin(
      earth["frame"].
      Translate(op(ListPadding(obj["coordinates"][1],3)))
      ));
  end if;

  out := plots:-display(
    plottools:-point(convert(O[1..3], list), symbol='solidbox', symbolsize = 20),
    linestyle = solid, color = ObjectColor(obj), scaling = constrained);

  PrintEndProc(procname);
  return out;
end proc: # PlotSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotDeformedSupport := proc(
  obj::{SUPPORT}, # Support to be plotted
  targets::{ # Support targets
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  {
    data::{list(`=`),set(`=`)} := [], # Substitutions
    scaling::{numeric}         := 1   # Scaling factor
  },
  $)::{procedure};

  description "Plot the deformed SUPPORT object <obj>.";

  local O, out;
  PrintStartProc(procname);
  #TODO: add compliant support deformation

  if (nops(obj["targets"]) > 1) then
    O := subs(op(data), Origin(
      GetObjByName(obj["targets"][2], targets)["frame"].
      Translate(op(ListPadding(obj["coordinates"][2],3)))
      ));
  else
    O := subs(op(data), Origin(
      earth["frame"].
      Translate(op(ListPadding(obj["coordinates"][1],3)))
      ));
  end if;

  out := plots:-display(
    plottools:-point(convert(O[1..3], list), symbol='solidbox', symbolsize = 20),
    linestyle = solid, color = ObjectColor(obj), parse("scaling") = constrained);

  PrintEndProc(procname);
  return out;
end proc: # PlotDeformedSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotStructure := proc(
  str::{STRUCTURE},                  # Structure to be plotted
  {
    data::{list(`=`),set(`=`)} := [] # Substitutions
  },
  $)::{list(procedure)};

  description "Plot the STRUCTURE object <str> given a list of substitutions <data>.";

  local out, disp, rb_joints, rb_loads, obj;
  PrintStartProc(procname);

  disp := []:
  for obj in str["objects"] do
    if IsBeam(obj) then
      disp := disp union [PlotBeam(obj, parse("data") = data)];
    elif IsRod(obj) then
      disp := disp union [PlotRod(obj, parse("data") = data)];
    elif IsSupport(obj) then
      disp := disp union [PlotSupport(obj, map(GetObjByName, obj["targets"], str["objects"]), parse("data") = data)];
    elif IsJoint(obj) then
      disp := disp union [PlotJoint(obj, map(GetObjByName, obj["targets"], str["objects"]), parse("data") = data)];
    elif IsRigidBody(obj) then
      GetObjsByType(['JOINT', 'SUPPORT'], str["objects"]);
      rb_joints := remove(x -> (not member(obj["name"], x["targets"])), %);
      GetObjsByType(['FORCE', 'MOMENT'], str["external_actions"]);
      rb_loads := remove(x -> obj["name"] <> x["target"], %);
      disp := disp union [PlotRigidBody(obj, rb_joints, rb_loads, parse("data") = data)];
    end if;
  end do;

  out := plots:-display(disp, axes = boxed, scaling = constrained, labels=['x', 'y', 'z']);

  PrintEndProc(procname);
  return out;
end proc: # PlotStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotDeformedStructure := proc(
  str::{STRUCTURE},                   # Structure to be plotted
  {
    data::{list(`=`),set(`=`)} := [], # Substitutions
    scaling::{numeric}         := 1   # Scaling factor
  },
  $)::{list(procedure)};

  description "Plot the deformed STRUCTURE object <str> given a list of "
              "substitutions <data> and a scaling factor <scaling>.";

  local out, disp, rb_joints, rb_loads, obj;
  PrintStartProc(procname);

  # Check if displacements and frame displacements are solved
  if (str["displacements_solved"] = false) or
      (str["frame_displacements_solved"] = false) then
    error("Displacements and frame displacements must be solved before plotting the deformed structure");
  end if;

  disp := []:
  for obj in str["objects"] do
    if IsBeam(obj) then
      disp := disp union [PlotDeformedBeam(obj, parse("data") = data, parse("scaling") = scaling)];
    elif IsRod(obj) then
      disp := disp union [PlotDeformedRod(obj, parse("data") = data, parse("scaling") = scaling)];
    elif IsSupport(obj) then
      disp := disp union [PlotDeformedSupport(obj, map(GetObjByName, obj["targets"], str["objects"]), parse("data") = data, parse("scaling") = scaling)];
    elif IsJoint(obj) then
      disp := disp union [PlotDeformedJoint(obj, map(GetObjByName, obj["targets"], str["objects"]), parse("data") = data, parse("scaling") = scaling)];
    elif IsRigidBody(obj) then
      GetObjsByType(['JOINT', 'SUPPORT'], str["objects"]);
      rb_joints := remove(x -> (not member(obj["name"], x["targets"])), %);
      GetObjsByType(['FORCE', 'MOMENT'], str["external_actions"]);
      rb_loads := remove(x -> obj["name"] <> x["target"], %);
      disp := disp union [PlotDeformedRigidBody(obj, rb_joints, rb_loads, parse("data") = data, parse("scaling") = scaling)];
    end if;
  end do;

  out := plots:-display(disp, axes = boxed, parse("scaling") = constrained, labels=['x', 'y', 'z']);

  PrintEndProc(procname);
  return out;
end proc: # PlotDeformedStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsInsideJoint := proc(
  obj::{JOINT},           # Joint object
  p::{POINT},             # Point to be checked
  tol::{numeric} := 1e-3, # Tolerance
  $)::{boolean};

  description "Check if the point <p> is inside the JOINT <obj>.";

  local O, out;
  PrintStartProc(procname);

  if not (nops(p) = 3) then
    error "The input point must be a list of 3 elements";
  end if;

  if (nops(obj["targets"]) > 1) then
    O := Origin(
      GetObjByName(obj["targets"][1], targets)["frame"].
      Translate(obj["coordinates"][1], 0, 0)
      );
  elif (not suppress_warnings) then
    WARNING("The support has no targets");
  end if;

  out := (norm(p - O) <= tol);

  PrintEndProc(procname);
  return out;
end proc: # IsInsideJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsInsideSupport := proc(
  obj::{SUPPORT},         # Support object
  p::{POINT},             # Point to be checked
  tol::{numeric} := 1e-3, # Tolerance
  $)::{boolean};

  description "Check if the point <p> is inside the SUPPORT <obj>.";

  local O, out;
  PrintStartProc(procname);

  if not (nops(p) = 3) then
    error "The input point must be a list of 3 elements";
  end if;

  if (nops(obj["targets"]) > 1) then
    O := Origin(
      GetObjByName(obj["targets"][2], targets)["frame"].
      Translate(obj["coordinates"][2], 0, 0)
      );
  elif (not suppress_warnings) then
    WARNING("The support has no targets");
  end if;

  out := (norm(p - O) <= tol);

  PrintEndProc(procname);
  return out;
end proc: # IsInsideSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsInsideRod := proc(
  obj::{ROD}, # Rod object
  p::{POINT}, # Point to be checked
  $)::{boolean};

  description "Check if the point <p> is inside the ROD <obj>.";

  local O, V, W, out;
  PrintStartProc(procname);

  if not (nops(p) = 3) then
    error "The input point must be a list of 3 elements";
  end if;

  O := Origin(obj["frame"]);
  V := obj["frame"].Translate(obj["length"], 0, 0) - O;
  W := p - O;

  out := (dot(W, V) >= 0) and (dot(W, V) <= dot(V, V));

  PrintEndProc(procname);
  return out;
end proc: # IsInsideRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsInsideBeam := proc(
  obj::{BEAM}, # Beam object
  p::{POINT},  # Point to be checked
  $)::{boolean};

  description "Check if the point <p> is inside the BEAM <obj>.";

  local O, V, W, out;
  PrintStartProc(procname);

  if not (nops(p) = 3) then
    error "The input point must be a list of 3 elements";
  end if;

  O := Origin(obj["frame"]);
  V := obj["frame"].Translate(obj["length"], 0, 0) - O;
  W := p - O;

  out := (dot(W, V) >= 0) and (dot(W, V) <= dot(V, V));

  PrintEndProc(procname);
  return out;
end proc: # IsInsideBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsInsideStructure := proc(
  obj::{STRUCTURE}, # Structure object
  p::{POINT},       # Point to be checked
  $)::{boolean};

  description "Check if the point <p> is inside the STRUCTURE <obj>.";

  local i, out;
  PrintStartProc(procname);

  out := false;
  for i in str["objects"] do
    if IsJoint(i) then
      out := out or IsInsideJoint(i, p);
    elif IsSupport(i) then
      out := out or IsInsideSupport(i, p);
    elif IsBeam(i) then
      out := out or IsInsideBeam(i, p);
    elif IsRod(i) then
      out := out or IsInsideRod(i, p);
    else
      error "Unknown object type";
    end if;
  end do;

  PrintEndProc(procname);
  return out;
end proc: # IsInsideStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module:

# That's all folks!
