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

TrussMe := module()

export  `union`,
        SetVerboseMode,
        PrintStartProc,
        PrintEndProc,
        IsEarth,
        SetGravity,
        GetGravity,
        Show,
        Rotate,
        Translate,
        Project,
        InverseFrame,
        IsFrame,
        Origin,
        Uvec,
        UvecX,
        UvecY,
        UvecZ,
        MakeMaterial,
        IsMaterial,
        MakeBeam,
        IsBeam,
        MakeRod,
        IsRod,
        MakeRigidBody,
        IsRigidBody,
        MakeJoint,
        IsJoint,
        MakeSupport,
        IsSupport,
        IsCompliantSupport,
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
        CleanJoint,
        CleanSupport,
        CleanRod,
        CleanBeam,
        CleanStructure,
        DrawStructureGraph;

global  ground;

local   ModuleLoad,
        ModuleUnload,
        earth,
        gravity,
        GetNames,
        ComputeDOF,
        NewtonEuler,
        HyperstaticSolver,
        IsostaticSolver,
        ComputeInternalActions,
        ComputePotentialEnergy,
        ComputeDisplacements,
        InternalActions,
        InitTrussMe,
        TypeRegister,
        Protect,
        Norm2,
        ComputeSpringDisplacement,
        ComputeSpringEnergy,
        ComputeSupportDisplacements,
        ComputeSupportInducedDisplacements,
        ObjectColor,
        PlotBeam,
        PlotRod,
        PlotJoint,
        PlotSupport,
        IsInsideJoint,
        IsInsideSupport,
        IsInsideBeam,
        IsInsideRod,
        IsInsideStructure,
        lib_base_path,
        verbose_mode,
        print_indent,
        print_increment,
        ListPadding,
        Beam_color,
        Rod_color,
        RigidBody_color,
        CompliantSupport_color,
        Support_color,
        CompliantJoint_color,
        Joint_color,
        Earth_color;

option  package,
        load   = ModuleLoad,
        unload = ModuleUnload;

description "A Maple Library for Truss Elements Structures";

#uses ; # add here the libraries used by the module

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   __  __           _       _      _                    _
#  |  \/  | ___   __| |_   _| | ___| |    ___   __ _  __| |
#  | |\/| |/ _ \ / _` | | | | |/ _ \ |   / _ \ / _` |/ _` |
#  | |  | | (_) | (_| | |_| | |  __/ |__| (_) | (_| | (_| |
#  |_|  |_|\___/ \__,_|\__,_|_|\___|_____\___/ \__,_|\__,_|

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ModuleLoad := proc()

  description "Module 'TrussMe' module load procedure";

  local i;

  # Display module init message
  printf(cat(
    "Module 'TrussMe' version beta-0.0, ",
    "Copyright (C) 2022-2023, M. Larcher & D. Stocco, ",
    "University of Trento 2022-2023"
    ));

  # Library path
  lib_base_path := null;
  for i in [libname] do
    if (StringTools[Search]("TrussMe", i) <> 0) then
      lib_base_path := i;
    end;
  end;
  if (lib_base_path = null) then
    error "cannot find 'TrussMe' library" ;
  end:

  # Register types
  TypeRegister();

  # Initialize the module variables
  InitTrussMe();

  # Protect Module Keywords
  Protect();

  NULL;
end proc: # ModuleLoad

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   __  __           _       _      _   _       _                 _
#  |  \/  | ___   __| |_   _| | ___| | | |_ __ | | ___   __ _  __| |
#  | |\/| |/ _ \ / _` | | | | |/ _ \ | | | '_ \| |/ _ \ / _` |/ _` |
#  | |  | | (_) | (_| | |_| | |  __/ |_| | | | | | (_) | (_| | (_| |
#  |_|  |_|\___/ \__,_|\__,_|_|\___|\___/|_| |_|_|\___/ \__,_|\__,_|

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ModuleUnload := proc()
  description "Module 'TrussMe' module unload procedure";
  printf("Unloading 'TrussMe'\n");
  verbose_mode := 1;

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

  description "Register 'TrussMe' module types";

  # Register types
  TypeTools[AddType](EARTH, IsEarth);
  TypeTools[AddType](FRAME, IsFrame);
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

  description "Initialize 'TrussMe' module internal variables";

  # Define Module Variables
  ground := <<1, 0, 0, 0>|
             <0, 1, 0, 0>|
             <0, 0, 1, 0>|
             <0, 0, 0, 1>>;

  gravity := [0, 0, 0];

  verbose_mode           := 1;
  print_indent           := 0;
  print_increment        := 4;
  Beam_color             := "SteelBlue";
  Rod_color              := "Niagara DarkOrchid";
  RigidBody_color        := "DarkKhaki";
  CompliantSupport_color := "DarkGreen";
  Support_color          := "DarkOrange";
  CompliantJoint_color   := "LightSalmon";
  Joint_color            := "MediumSeaGreen";
  Earth_color            := "Firebrick";

  earth := table({
    parse("type")             = EARTH,
    parse("name")             = "earth",
    parse("length")           = 0,
    parse("frame")            = ground,
    parse("admissible_loads") = [1, 1, 1, 1, 1, 1]
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
  B::{list, set},  # Object B to be united
  {
    verbose::boolean := false # Verbose mode
  }, $) option overload;

  description "Extension of union operator to list objects <A> and <B>";

  local out;
  PrintStartProc(procname);

  if type(A, 'set') and type(B, 'set') then
    out := {op(A), op(B)};
  else
    out := [op(A), op(B)];
  end if:

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

SetVerboseMode := proc(
  mode::integer, # Procedure name
  $)::nothing;

  description "Set verbose mode <mode> between 0, 1, 2";

  # Show start message
  if (mode < 0) or
     (mode > 2) then
    error "invalid verbose mode detected";
  else
    verbose_mode := mode;
  end if;
  NULL;
end proc: # PrintStartProc

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PrintStartProc := proc(
  proc_name::procedure, # Procedure name
  $)::nothing;

  description "Print the start message of a procedure with name <proc_name>";

  # Increase printf indentation
  print_indent := print_indent + print_increment;

  # Show start message
  if (verbose_mode > 1) then
    printf("%*sStart '%s' procedure...\n", print_indent, "", proc_name);
  end if:
end proc: # PrintStartProc


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PrintEndProc := proc(
  proc_name::procedure, # Procedure name
  $)::nothing;

  description "Print the end message of a procedure with name <proc_name>";

  # Show end message
  if (verbose_mode > 1) then
    printf("%*sEnd   '%s' procedure\n", print_indent, "", proc_name);
  end if:

  # Increase printf indentation
  print_indent := print_indent - print_increment;
end proc: # PrintEndProc

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsEarth := proc(
  obj::anything, # Object to be tested
  $)::boolean;

  description "Test if an object <obj> is the EARTH object";

  local out;
  PrintStartProc(procname);

  if type(obj, table) and
     (obj[parse("type")] = EARTH) and
     (obj[parse("length")] = 0) and
     (obj[parse("frame")] = ground) and
     (obj[parse("admissible_loads")] = [1, 1, 1, 1, 1, 1]) then
    out := true;
  else
    out := false;
  end if:

  PrintEndProc(procname);
  return out;
end proc: # IsEarth

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SetGravity := proc(
  obj::list, # Gravity vector
  $)::nothing;

  description "Set gravity vector with [x, y, z] components of <obj>";

  PrintStartProc(procname);
  # Set gravity local variable
  if (nops(obj) = 3) then
    gravity := obj;
  else
    error "invalid gravity vector detected";
  end if:
  PrintEndProc(procname);
  NULL;
end proc: # SetGravity

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

GetGravity := proc(
  $)::list;

  description "Get gravity vector";

  PrintStartProc(procname);
  PrintEndProc(procname);
  return gravity;
end proc: # GetGravity

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Norm2 := proc(
  obj::{list, vector}, # Vector for which the norm is computed
  $)::algebraic;

  description "Compute the norm of a vector <obj>";

  local out;
  PrintStartProc(procname);

  out := sqrt(add(x, x in obj^~2));

  PrintEndProc(procname);
  return out;
end proc: # Norm2

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ListPadding := proc(
  list::{list, algebraic}, # List to be padded
  n::integer,              # Number of elements of the final list
  value::algebraic := 0,   # Value to be used for padding
$)::list;

  description "Pad a list <list> with <value> to have <n> elements";

  local i, out;
  PrintStartProc(procname);

  if type(list, algebraic) then
    out := [list];
  else
    out := list;
  end if:

  if (nops(out) < n) then
    out := out union [seq(value, i = (1..n-nops(out)))];
  end if:

  PrintEndProc(procname);
  return out;
end proc: # ListPadding

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Show := proc(
  tab::table, # Table to be shown
  $)::nothing;

  description "Show the content of a table <tab>";

  PrintStartProc(procname);
  print(tab = tab[parse("type")](op(op(tab))));
  PrintEndProc(procname);
end proc: # Show

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

GetNames := proc(
  objs::{ # Structural elements
    list({MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH}),
    set( {MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH})
  }, $)::{list({string}), set({string})};

  description "Get names of a list/set of objects <objs>";

  local out;
  PrintStartProc(procname);

  if type(objs, 'set') then
    out := {seq(objs[i][parse("name")], i = 1..nops(objs))};
  else
    out := [seq(objs[i][parse("name")], i = 1..nops(objs))];
  end if:

  PrintEndProc(procname);
  return out;
end proc: # GetNames

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

InverseFrame := proc(
  RF::FRAME, # Reference frame (affine transformation) to be inverted
  $)::FRAME;

  description "Inverse transformation matrix of an affine transformation <RF>";

  return LinearAlgebra:-Transpose(RF);
end proc: # InverseFrame

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsFrame := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the input object <obj> is a FRAME object";

  local out;
  PrintStartProc(procname);

  if (type(obj, 'Matrix')) and
     (LinearAlgebra:-RowDimension(obj) = 4) and
     (LinearAlgebra:-ColumnDimension(obj) = 4) and
     (simplify(combine(LinearAlgebra:-Determinant(obj) = 1))) then
    out := true;
  else
    out := false;
  end if;

  PrintEndProc(procname);
  return out;
end proc: # IsFrame

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Rotate := proc(
  axis::symbol,     # Rotation axis
  angle::algebraic, # Rotation angle (rad)
  $)::FRAME;

  description "Transformation matrix corresponding to the rotation <angle> "
    "around the given <axis>";

  local out;
  PrintStartProc(procname);

  if (axis = 'X') then
    out := <<1, 0,           0,          0>|
            <0, cos(angle),  sin(angle), 0>|
            <0, -sin(angle), cos(angle), 0>|
            <0, 0,           0,          1>>;
  elif (axis = 'Y') then
    out := <<cos(angle), 0, -sin(angle), 0>|
            <0,          1, 0,           0>|
            <sin(angle), 0, cos(angle),  0>|
            <0,          0, 0,           1>>;
  elif (axis = 'Z') then
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
  x::algebraic, # X-axis translation component
  y::algebraic, # Y-axis translation component
  z::algebraic, # Z-axis translation component
  $)::FRAME;

  description "Transformation matrix corresponding to the translation <x,y,z>";

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
  RF::FRAME, # Reference frame
  $)::vector;

  description "Extract the origin of the reference frame <RF>";

  local out;
  PrintStartProc(procname);

  out := <RF[1,4], RF[2,4], RF[3,4], 1>;

  PrintEndProc(procname);
  return out;
end proc: # Translate

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Uvec := proc(
  axis::symbol,        # Axis of the unit vector
  RF::FRAME := ground, # Reference frame
  $)::vector;

  description "Extract the unit vector of the reference frame <RF> along the "
    "given <axis>";

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
  end if:

  PrintEndProc(procname);
  return out;
end proc: # Uvec

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

UvecX := proc(
  RF::FRAME := ground, # Reference frame
  $)::vector;

  description "Extract the x-axis unit vector of the reference frame <RF>";

  local out;
  PrintStartProc(procname);

  out := <RF[1,1], RF[2,1], RF[3,1], 0>;

  PrintEndProc(procname);
  return out;
end proc: # UvecX

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

UvecY := proc(
  RF::FRAME := ground, # Reference frame
  $)::vector;

  description "Extract the y-axis unit vector of the reference frame <RF>";

  local out;
  PrintStartProc(procname);

  out := <RF[1,2], RF[2,2], RF[3,2], 0>;

  PrintEndProc(procname);
  return out;
end proc: # UvecY

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

UvecZ := proc(
  RF::FRAME := ground, # Reference frame
  $)::vector;

  description "Extract the z-axis unit vector of the reference frame <RF>";

  local out;
  PrintStartProc(procname);

  out := <RF[1,3], RF[2,3], RF[3,3], 0>;

  PrintEndProc(procname);
  return out;
end proc: # UvecZ

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Project := proc(
  x::{list, vector}, # Object to be projected
  RF_ini::FRAME,     # Reference frame from which the object is expressed
  RF_end::FRAME,     # Reference frame to which the object will be expressed
  $)::{list, vector};

  description "Project <x,y,z>, or vector <x,y,z,0>, or point <x,y,z,1> from "
    "reference frame <RF_ini> to reference frame <RF_end>";

  local x_tmp, out;

  PrintStartProc(procname);

  # Pad input vector with 0 if its length is 3
  if (nops(x) = 3) then
    x_tmp := <x[1], x[2], x[3], 0>;
  elif (nops(x) = 4) then
    x_tmp := <x[1], x[2], x[3], x[4]>;
  else
    error "invalid input vector/point <x> detected";
  end if;

  LinearAlgebra[MatrixInverse](RF_end).RF_ini.x_tmp;
  out := simplify([seq(%[i], i = 1..nops(x))]);

  PrintEndProc(procname);
  return out;
end proc: # Project

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeMaterial := proc({
    name::string               := "steel",   # Name of the material
    elastic_modulus::algebraic := 210.0E+09, # Elastic modulus (Pa)
    poisson_ratio::algebraic   := 0.3,       # Poisson ratio (-)
    shear_modulus::algebraic   := elastic_modulus/(2*(1+poisson_ratio)),
                                             # Shear modulus (Pa)
    density::algebraic         := 7.4E+03    # Density (kg/m^3)
  }, $)::MATERIAL;

  description "Define a MATERIAL object with inputs: name of the material, "
    "elastic modulus <elastic_modulus> (default = 210.0E9 Pa), Poisson ratio "
    "<poisson_ratio> (default = 0.3), shear modulus <shear_modulus> (default "
    "= E/(2*(1+nu))), density <density> (default = 7.4E3 kg/m^3)";

  local out;
  PrintStartProc(procname);

  out := table({
    parse("type")            = MATERIAL,
    parse("name")            = name,
    parse("elastic_modulus") = elastic_modulus,
    parse("poisson_ratio")   = poisson_ratio,
    parse("shear_modulus")   = shear_modulus,
    parse("density")         = density
    });

  PrintEndProc(procname);
  return op(out);
end proc: # DefineMaterial

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsMaterial := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the input object <obj> is a MATERIAL object";

  local out;
  PrintStartProc(procname);

  if (obj[parse("type")] = MATERIAL) and
     type(obj, table) and
     type(obj[parse("name")], string) and
     type(obj[parse("elastic_modulus")], algebraic) and
     type(obj[parse("poisson_ratio")], algebraic) and
     type(obj[parse("shear_modulus")], algebraic) and
     type(obj[parse("density")], algebraic) then
    out := true;
  else
    out := false;
  end if;

  PrintEndProc(procname);
  return out;
end proc: # IsMaterial

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeForce := proc(
  components::list,                                    # Force components in RF
  coords::{algebraic, list(algebraic)},                # Application coordinates in object frame
  obj::{BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH}, # Target object
  RF::FRAME := ground,                                 # Reference frame
  $)::FORCE;

  description "Define a FORCE object with inputs: force components <components>, "
    "force application axial coordinate <coords>, target object <obj>, and "
    "optional reference frame <RF> in which the force is defined (default = ground)";

  local proj_components, out;

  PrintStartProc(procname);

  proj_components := Project(components, RF, obj[parse("frame")]);
  # Check input arguments
  if IsRod(obj) then
    if (proj_components[2] <> 0) or (proj_components[3] <> 0) then
      error "only axial forces are accepted in ROD objects";
    end if;
  elif IsSupport(obj) or IsJoint(obj) then
    if (evalf(coords) <> 0) then
      error "only null axial coordinate is accepted for SUPPORT and JOINT "
        "objects";
    end if;
  end if;

  out := table({
    parse("type")       = FORCE,
    parse("components") = proj_components,
    parse("coordinate") = ListPadding(coords, 3),
    parse("target")     = obj[parse("name")]
    });

  PrintEndProc(procname);
  return op(out);
end proc: # MakeForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsForce := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a FORCE object";

  local out;
  PrintStartProc(procname);

  if (obj[parse("type")] = FORCE) and
     type(obj, table) and
     type(obj[parse("components")], list) and
     type(obj[parse("coordinate")], {algebraic, list(algebraic)}) and
     type(obj[parse("target")], string) then
    out := true;
  else
    out := false;
  end if;

  PrintEndProc(procname);
  return out;
end proc: # IsForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeMoment := proc(
  components::list,                               # Moment components
  coords::{algebraic, list(algebraic)},           # Application coordinates in object frame
  obj::{BEAM, RIGID_BODY, SUPPORT, JOINT, EARTH}, # Target object
  RF::FRAME := ground,                            # Reference frame
  $)::MOMENT;

  description "Define a MOMENT object with inputs: moment components <components>, "
    "moment application axial coordinate <coords>, target object <obj>, and "
    "optional reference frame <RF> in which the moment is  defined (default "
    "= ground)";

  local proj_components, out;
  PrintStartProc(procname);

  proj_components := Project(components, RF, obj[parse("frame")]);
  # Check input arguments
  if IsSupport(obj) or IsJoint(obj) then
    if (evalf(coords) <> 0) then
      error "only null axial coordinate is accepted for SUPPORT and JOINT "
        "objects";
    end if;
  end if;

  out := table({
    parse("type")       = MOMENT,
    parse("components") = proj_components,
    parse("coordinate") = ListPadding(coords, 3),
    parse("target")     = obj[parse("name")]
    });

  PrintEndProc(procname);
  return op(out);
end proc: # MakeMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsMoment := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a MOMENT object";

  local out;
  PrintStartProc(procname);

  if (obj[parse("type")] = MOMENT) and
     type(obj, table) and
     type(obj[parse("components")], list) and
     type(obj[parse("coordinate")], {algebraic, list(algebraic)}) and
     type(obj[parse("target")], string) then
    out := true;
  else
    out := false;
  end if;

  PrintEndProc(procname);
  return out;
end proc: # IsMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeQForce := proc(
  components::{procedure,list(algebraic)}, # Distributed load components
  obj::{BEAM, ROD},                        # Target object
  RF::FRAME := ground,                     # Reference frame
  {
    ell_min::algebraic := 0,                   # Initial axial coordinate
    ell_max::algebraic := obj[parse("length")] # Final axial coordinate
  }, $)::QFORCE;

  description "Define a QFORCE object with inputs: distributed load target "
    "object components <components>, target object <obj>, optional reference "
    "frame <RF> in which the load components are defined (default = ground), "
    "and optional initial <ell_min> and final <ell_max> application points "
    "(axial coordinates)";

  local proj_components, x, out;
  PrintStartProc(procname);

  if type(components, procedure) then
    proj_components := unapply(Project(components(x), RF, obj[parse("frame")]), x);
  else
    proj_components := (x) -> piecewise((ell_min <= x) and (x <= ell_max), Project(components, RF, obj[parse("frame")]), 0);
  end if;

  if IsRod(obj) then
    if (proj_components(x)[2] <> 0) or (proj_components(x)[3] <> 0) then
      error "only axial loads are accepted in ROD objects"
    end if;
  end if;

  out := table({
    parse("type")        = QFORCE,
    parse("components")  = proj_components,
    parse("target")      = obj[parse("name")]
    });

  PrintEndProc(procname);
  return op(out);
end proc: # MakeQForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsQForce := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a QFORCE object";

  local out;
  PrintStartProc(procname);

  if (obj[parse("type")] = QFORCE) and
     type(obj, table) and
     type(obj[parse("components")], procedure) and
     type(obj[parse("target")], string) then
    out := true;
  else
    out := false;
  end if;

  PrintEndProc(procname);
  return out;
end proc: # IsQForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeQMoment := proc(
  components::{procedure,list(algebraic)}, # Distributed load components
  obj::BEAM,                               # Target object
  RF::FRAME := ground,                     # Reference frame in which the moment is defined
  {
    ell_min::algebraic := 0,                   # Initial application point (axial coordinate)
    ell_max::algebraic := obj[parse("length")] # Final application point (axial coordinate)
  }, $)::QMOMENT;

  description "Define a QMOMENT object with inputs: distributed torque target "
    "object components <components>, target object <obj>, optional reference "
    "frame <RF> in which the load components are defined (default = ground), "
    "and optional initial <ell_min> and final <ell_max> application points "
    "(axial coordinates)";

  local proj_components, x, out;
  PrintStartProc(procname);

  if type(components, procedure) then
    proj_components := unapply(Project(components(x), RF, obj[parse("frame")]),x);
  else
    proj_components := (x) -> piecewise((ell_min <= x) and (x <= ell_max), Project(components, RF, obj[parse("frame")]), 0);
  end if;

  out := table({
    parse("type")        = QMOMENT,
    parse("components")  = proj_components,
    parse("target")      = obj[parse("name")]
    });

  PrintEndProc(procname);
  return op(out);
end proc: # MakeQMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsQMoment := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a QMOMENT object";

  local out;
  PrintStartProc(procname);

  if (obj[parse("type")] = QMOMENT) and
     type(obj, table) and
     type(obj[parse("components")], procedure) and
     type(obj[parse("target")], string) then
    out := true;
  else
    out := false;
  end if;

  PrintEndProc(procname);
  return out;
end proc: # IsQMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeSupport := proc(
  name::string,                        # Support name
  constrained_dof::list,               # Constrained degree of freedom
  objs::list({BEAM, ROD, RIGID_BODY}), # Target objects
  coords::list,                        # Support locations
  RF::FRAME := ground,                 # Reference frame of the support
  {
    stiffness::{procedure,list(algebraic)} := [ # Stiffness components (default = infinite)
      infinity, infinity, infinity,
      infinity, infinity, infinity
    ] *~ constrained_dof
  }, $)::SUPPORT;

  description "Make a SUPPORT object with inputs: support name <name>, constrained "
    "degrees of freedom <constrained_dof>, target objects <objs>, support locations "
    "<coords>, and optional reference frame <RF> in which the support is defined "
    "(default = ground)";

  local S, J_tmp, i, j, sr_F_names, sr_F_values_tmp, sr_M_names, sr_M_values_tmp,
    S_stiffness, x;
  PrintStartProc(procname);

  for i from 1 to nops(objs) do
    if IsRod(objs[i]) and (coords[i] <> 0) and (coords[i] <> objs[i][parse("length")]) then
      error "SUPPORT objects can only be applied at extremes of ROD objects"
    end if;
    if IsRod(objs[i]) and (constrained_dof[4..6] <> [0, 0, 0]) then
      error "ROD objects supports can only have translational constraints"
    end if;
  end do;

  if type(stiffness, procedure) then
    S_stiffness := unapply(stiffness(x) *~ constrained_dof, x);
    # Check for non zero stiffness on constrained dof
    if has(stiffness(x)[remove(x -> x=0, ([seq(i, i = 1..6)]) *~ constrained_dof)], 0) then
      error "stiffness corresponding to constrained degrees of freedom cannot be zero";
    end if;
    # Check for zero stiffness on unconstrained dof
    if (S_stiffness(x) <> stiffness(x)) then
      WARNING("stiffness components not corresponding to constrained_dof are ignored");
    end if;
  else
    # Check for non zero stiffness on constrained dof
    if has(stiffness[remove(x -> x=0, ([seq(i, i = 1..6)]) *~ constrained_dof)], 0) then
      error "stiffness corresponding to constrained degrees of freedom cannot be zero";
    end if;
    # Check for zero stiffness on unconstrained dof
    S_stiffness := (x) -> stiffness *~ constrained_dof;
    if (S_stiffness(x) <> stiffness) then
      WARNING("stiffness components not corresponding to constrained_dof are ignored");
    end if;
  end if;

  S := table({
    parse("type")                     = SUPPORT,
    parse("constrained_dof")          = constrained_dof,
    parse("coordinates")              = [[0,0,0], op(map(ListPadding, coords, 3))],
    parse("name")                     = name,
    parse("frame")                    = RF,
    parse("targets")                  = [earth[parse("name")], op(GetNames(objs))],
    parse("variables")                = [],
    parse("forces")                   = [],
    parse("moments")                  = [],
    parse("constraint_loads")         = [],
    parse("support_reactions")        = [], # Expressed in support reference frame
    parse("stiffness")                = S_stiffness,
    parse("displacements")            = []
    });

  # Build the temporary joint
  J_tmp := MakeJoint(name, constrained_dof, [earth, op(objs)], S[parse("coordinates")], RF);

  S[parse("variables")]                := J_tmp[parse("variables")];
  S[parse("forces")]                   := J_tmp[parse("forces")];
  S[parse("moments")]                  := J_tmp[parse("moments")];
  S[parse("constraint_loads")]         := J_tmp[parse("constraint_loads")];

  # Retrieve support force reactions
  sr_F_names := [FX, FY, FZ];
  for i from 1 to nops(S[parse("forces")]) do
    if (S[parse("forces")][i][parse("target")] = earth[parse("name")]) then
      # Project forces in the support reference frame
      sr_F_values_tmp := Project(S[parse("forces")][i][parse("components")], ground, S[parse("frame")]);
      for j from 1 to 3 do
        if (sr_F_values_tmp[j] <> 0) then
          S[parse("support_reactions")] := [
            op(S[parse("support_reactions")]),
            sr_F_names[j] = -sr_F_values_tmp[j]
            ];
        end if;
      end do;
      break;
    end if;
  end do;

  # Retrieve support moments reactions
  sr_M_names := [MX, MY, MZ];
  for i from 1 to nops(S[parse("moments")]) do
    if (S[parse("moments")][i][parse("target")] = earth[parse("name")]) then
      # Project moments in the support reference frame
      sr_M_values_tmp := Project(S[parse("moments")][i][parse("components")], ground, S[parse("frame")]);
      for j from 1 to 3 do
        if (sr_M_values_tmp[j] <> 0) then
          S[parse("support_reactions")] := [
            op(S[parse("support_reactions")]),
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
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a SUPPORT object";

  local out;
  PrintStartProc(procname);

  if (obj[parse("type")] = SUPPORT) and
     type(obj[parse("constrained_dof")], list) and
     type(obj[parse("coordinates")], list) and
     type(obj[parse("name")], string) and
     type(obj[parse("frame")], FRAME) and
     type(obj[parse("targets")], list({string})) and
     type(obj[parse("variables")], list) and
     type(obj[parse("forces")], list) and
     type(obj[parse("moments")], list) and
     type(obj[parse("constraint_loads")], list) and
     type(obj[parse("support_reactions")], list) and
     type(obj[parse("stiffness")], procedure) and
     type(obj[parse("displacements")], list) then
    out := true;
  else
    out := false;
  end if;

  PrintEndProc(procname);
  return out;
end proc: # IsSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsCompliantSupport := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a SUPPORT object with compliant "
    "constraints";

  local out, i;
  PrintStartProc(procname);

  out := false;
  if IsSupport(obj) then
    for i from 1 to nops(obj[parse("stiffness")]) do
      if (obj[parse("stiffness")](x)[i] <> infinity) and
         (obj[parse("constrained_dof")][i] = 1) then
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
  obj::SUPPORT, # Support to be cleaned
  $)::nothing;

  description "Clean SUPPORT object <obj> internal variables";

  PrintStartProc(procname);
  obj[parse("constraint_loads")]         := [];
  obj[parse("support_reactions")]        := [];
  obj[parse("displacements")]            := [];
  PrintEndProc(procname);
end proc: # CleanSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeJoint := proc(
  name::string,                                               # Joint name
  constrained_dof::list,                                      # Constrained degree of freedom
  objs::list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH}), # Target objects
  coords::list,                                               # Joint locations
  RF::FRAME := ground,                                        # Reference frame
  $)::JOINT;

  description "Make a JOINT object with inputs: joint name <name>, constrained "
    "degrees of freedom <constrained_dof>, target objects <objs>, joint locations "
    "<coords>, and optional reference frame <RF> in which the joint is defined "
    "(default = ground)";

  local J, i, jf_comp, jm_comp, jf_comp_obj, jm_comp_obj, jm_indets, jf_indets,
    constraint, P_tmp;
  PrintStartProc(procname);

  for i from 1 to nops(objs) do
    if IsRod(objs[i]) and
        (ListPadding(coords[i], 3) <> [0, 0, 0]) and
        (ListPadding(coords[i], 3) <> [objs[i][parse("length")], 0, 0]) then
      error "JOINT objects can only be applied at extremes of ROD objects";
    end if;
    if IsRod(objs[i]) and (constrained_dof[4..6] <> [0, 0, 0]) then
      error "ROD objects supports can only have translational constraints";
    end if;
  end do;

  J := table({
    parse("type")                     = JOINT,
    parse("constrained_dof")          = constrained_dof,
    parse("coordinates")              = coords,
    parse("name")                     = name,
    parse("frame")                    = RF,
    parse("targets")                  = GetNames(objs),
    parse("variables")                = [],
    parse("forces")                   = [],
    parse("moments")                  = [],
    parse("constraint_loads")         = []
    });

  # Check if joint position on each object is the same
  # FIXME: does not work for mixed numerical and symbolic coordinates
  #if nops(ells) > 1 then
  #  P_tmp := objs[1][parse("frame")].Translate(ells[1], 0, 0);
  #  for i from 2 to nops(ells) do
  #    if not (norm(P_tmp - objs[i][parse("frame")].Translate(ells[i], 0, 0)) = 0) then
  #      error "Joint locations are not the same on all objects";
  #    end if;
  #  end do;
  #end if;

  # Add all the bodies forces
  for i from 1 to nops(objs) do
    # Create force compatible with the joint constrained dof
    jf_comp := convert(<
      JFx_||(J[parse("name")])||_||(objs[i][parse("name")]),
      JFy_||(J[parse("name")])||_||(objs[i][parse("name")]),
      JFz_||(J[parse("name")])||_||(objs[i][parse("name")])
      > *~ <op(constrained_dof[1..3])>,
      list);
    # Project the components into object frame and extract admissible loads
    jf_comp_obj := convert(
      Project(jf_comp, RF, objs[i][parse("frame")])
      .~ <op(objs[i][parse("admissible_loads")][1..3])>,
      list);
    # Extract the survived components
    jf_indets := indets(jf_comp_obj);
    # Use the non admissible loads to build the loads constraint
    constraint := convert(
      Project(jf_comp, RF, objs[i][parse("frame")])
      .~ <op((-1*objs[i][parse("admissible_loads")][1..3]) +~ 1)>,
      list);
    # Remove the null equations
    constraint := remove(x -> x = 0, constraint);
    # Remove useless constraints
    constraint := select(x -> has(x, jf_indets), constraint);
    # Update the joint constraint loads
    J[parse("constraint_loads")] := J[parse("constraint_loads")] union constraint;
    # Check if there are reactions
    if (jf_comp_obj <> [0, 0, 0]) then
      # Create the reaction force between joint and obj
      JF_||(name)||_||(objs[i][parse("name")]) := MakeForce(
        jf_comp_obj, coords[i], objs[i], objs[i][parse("frame")]
        );
      JF_||(objs[i][parse("name")])||_||(name) := MakeForce(
        -jf_comp_obj, 0, J, objs[i][parse("frame")]
        );
      # Update the output joint
      J[parse("variables")] := J[parse("variables")] union jf_indets;
      J[parse("forces")] := [
        op(J[parse("forces")]),
        JF_||(name)||_||(objs[i][parse("name")]),
        JF_||(objs[i][parse("name")])||_||(name)
        ];
    end if;
  end do;

  # Add all the bodies moments
  for i from 1 to nops(objs) do
    # Create moment compatible with joint constrained dof
    jm_comp := convert(<
      JMx_||(J[parse("name")])||_||(objs[i][parse("name")]),
      JMy_||(J[parse("name")])||_||(objs[i][parse("name")]),
      JMz_||(J[parse("name")])||_||(objs[i][parse("name")])
      > *~ <op(constrained_dof[4..6])>,
      list);
    # Project the components into object frame and extract the admissible loads
    jm_comp_obj := convert(
      Project(jm_comp, RF, objs[i][parse("frame")])
      .~ <op(objs[i][parse("admissible_loads")][4..6])>,
      list);
    # Use the non admissible loads to build the loads constraint
    constraint := convert(
      Project(jm_comp, RF, objs[i][parse("frame")])
      .~ <op((-1*objs[i][parse("admissible_loads")][4..6]) +~ 1)>,
      list);
    constraint := remove(x -> x = 0, constraint);
    J[parse("constraint_loads")] := [
      op(J[parse("constraint_loads")]),
      op(constraint)
      ];
    # Extract the survived components
    jm_indets := indets(jm_comp);
    # Check if there are reactions
    if not ((add(jm_comp_obj) = 0) or (add(jm_comp_obj) = 0.)) then
      # Create the reaction force between joint and obj
      JM_||(name)||_||(objs[i][parse("name")]) := MakeMoment(
        jm_comp_obj, coords[i], objs[i], objs[i][parse("frame")]
        );
      JM_||(objs[i][parse("name")])||_||(name) := MakeMoment(
        -jm_comp_obj, 0, J, objs[i][parse("frame")]
        );
      # Update the output joint
      J[parse("variables")] := [
        op(J[parse("variables")]),
        op(jm_indets)
        ];
      J[parse("moments")] := [
        op(J[parse("moments")]),
        JM_||(name)||_||(objs[i][parse("name")]),
        JM_||(objs[i][parse("name")])||_||(name)
        ];
    end if;
  end do;

  PrintEndProc(procname);
  return op(J);
end proc: # MakeJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsJoint := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the objecy <obj> is a JOINT object";

  local out;
  PrintStartProc(procname);

  if (obj[parse("type")] = JOINT) and
     type(obj, table) and
     type(obj[parse("constrained_dof")], list) and
     type(obj[parse("coordinates")], list) and
     type(obj[parse("name")], string) and
     type(obj[parse("frame")], FRAME) and
     type(obj[parse("targets")], list({string})) and
     type(obj[parse("variables")], list) and
     type(obj[parse("forces")], list) and
     type(obj[parse("moments")], list) and
     type(obj[parse("constraint_loads")], list) then
    out := true;
  else
    out := false;
  end if;

  PrintEndProc(procname);
  return out;
end proc: # IsJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanJoint := proc(
  obj::JOINT, # Object to be cleaned
  $)::nothing;

  description "Clean JOINT object <obj> internal variables";

  PrintStartProc(procname);
  obj[parse("constraint_loads")] := [];
  PrintEndProc(procname);
end proc: # CleanJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeRod := proc(
  name::string,        # Object name
  ell::algebraic,      # Length (m)
  RF::FRAME := ground, # Reference frame
  {
    area::{algebraic, procedure} := infinity, # Section area (m^2)
    material::MATERIAL           := NULL      # Material
  }, $)::ROD;

  description "Create a ROD object with inputs: object name <name>, reference "
  "length <ell>, optional reference frame <RF> in which the rod is defined, and "
  "optional section area <area> and material type <material>";

  local area_proc, out;
  PrintStartProc(procname);

  if type(area, procedure) then
    area_proc := area;
  else
    area_proc := (x) -> area;
  end if;

  out := table({
    parse("type")             = ROD,
    parse("name")             = name,
    parse("length")           = ell,
    parse("area")             = area_proc,
    parse("material")         = material,
    parse("frame")            = RF,
    parse("admissible_loads") = [1, 0, 0, 0, 0, 0],
    parse("internal_actions") = [],
    parse("displacements")    = []
    });

  PrintEndProc(procname);
  return op(out);
end proc: # MakeRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsRod := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a ROD object";

  local out;
  PrintStartProc(procname);

  if (obj[parse("type")] = ROD) and
     type(obj, table) and
     type(obj[parse("name")], string) and
     type(obj[parse("length")], algebraic) and
     type(obj[parse("area")], procedure) and
     type(obj[parse("material")], MATERIAL) and
     type(obj[parse("frame")], FRAME) and
     type(obj[parse("admissible_loads")], list) and
     type(obj[parse("internal_actions")], list) and
     type(obj[parse("displacements")], list) then
    out := true;
  else
    out := false;
  end if;

  PrintEndProc(procname);
  return out;
end proc: # IsRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanRod := proc(
  obj::ROD, # Object to be cleaned
  $)::nothing;

  description "Clean ROD object <obj> internal variables";

  PrintStartProc(procname);
  obj[parse("internal_actions")] := [];
  obj[parse("displacements")]    := [];
  PrintEndProc(procname);
end proc: # CleanRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeBeam := proc(
  name::string,        # Object name
  ell::algebraic,      # Length (m)
  RF::FRAME := ground, # Reference frame
  {
    area::{algebraic, procedure}                   := infinity,   # Section area (m^2)
    timo_shear_coeff::{list(algebraic), procedure} := [5/6, 5/6], # Timoshenko shear coefficient
    material::MATERIAL                             := NULL,       # Material object
    I_xx::{algebraic, procedure}                   := infinity,   # Section x-axis inertia (m^4)
    I_yy::{algebraic, procedure}                   := infinity,   # Section y-axis inertia (m^4)
    I_zz::{algebraic, procedure}                   := infinity    # Section z-axis inertia (m^4)
  }, $)::BEAM;

  description "Create a BEAM object with inputs: object name <name>, reference "
    "length <ell>, optional reference frame <RF> in which the rod is defined, and "
    "optional section area <area> and material type <material> and inertias on "
    "x- <I_xx>, y- <I_yy>, and z-axis <I_zz>";

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
    parse("type")             = BEAM,
    parse("name")             = name,
    parse("length")           = ell,
    parse("area")             = area_proc,
    parse("timo_shear_coeff") = timo_shear_coeff_proc,
    parse("material")         = material,
    parse("inertias")         = [I_xx_proc, I_yy_proc, I_zz_proc],
    parse("frame")            = RF,
    parse("admissible_loads") = [1, 1, 1, 1, 1, 1],
    parse("internal_actions") = [],
    parse("displacements")    = []
    });

  PrintEndProc(procname);
  return op(out);
end proc: # MakeBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsBeam := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a BEAM object";

  local out;
  PrintStartProc(procname);

  if (obj[parse("type")] = BEAM) and
     type(obj, table) and
     type(obj[parse("name")], string) and
     type(obj[parse("length")], algebraic) and
     type(obj[parse("area")], procedure) and
     type(obj[parse("timo_shear_coeff")], procedure) and
     type(obj[parse("material")], MATERIAL) and
     type(obj[parse("inertias")], list(procedure)) and
     type(obj[parse("frame")], FRAME) and
     type(obj[parse("admissible_loads")], list) and
     type(obj[parse("internal_actions")], list) and
     type(obj[parse("displacements")], list) then
    out := true;
  else
    out := false;
  end if;

  PrintEndProc(procname);
  return out;
end proc: # IsBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanBeam := proc(
  obj::BEAM, # Object to be cleaned
  $)::nothing;

  description "Clean BEAM object <obj> internal variables";

  PrintStartProc(procname);
  obj[parse("internal_actions")] := [];
  obj[parse("displacements")]    := [];
  PrintEndProc(procname);
end proc: # CleanBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeRigidBody := proc(
  name::string,        # Object name
  RF::FRAME := ground, # Reference frame
  {
    COM::list(algebraic) := Origin(RF), # COM position in RF (default: Origin(RF))
    mass::algebraic := 0                # Mass (kg)
  },
  $)::RIGID_BODY;

  description "Create a RIGID_BODY object with inputs: object name <name>, "
    "reference frame <RF> in which the rigid body is defined, and optional "
    "center of mass position <COM> and mass <mass>";

  local out;
  PrintStartProc(procname);

  out := table({
    parse("type") = RIGID_BODY,
    parse("name") = name,
    parse("frame") = RF,
    parse("COM") = COM,
    parse("mass") = mass,
    parse("admissible_loads") = [1, 1, 1, 1, 1, 1]
    });

  PrintEndProc(procname);
  return op(out);
end proc: # MakeRigidBody

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsRigidBody := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a RIGID_BODY object";

  local out;
  PrintStartProc(procname);

  if (obj[parse("type")] = RIGID_BODY) and
     type(obj, table) and
     type(obj[parse("name")], string) and
     type(obj[parse("frame")], FRAME) and
     type(obj[parse("COM")], list(algebraic)) and
     type(obj[parse("mass")], algebraic) and
     type(obj[parse("admissible_loads")], list) then
    out := true;
  else
    out := false;
  end if;

  PrintEndProc(procname);
  return out;
end proc: # IsRigidBody

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeSpringDisplacement := proc(
  spring_load::algebraic,      # Load on the spring
  spring_stiffness::procedure, # Spring stiffness
  $)::algebraic;

  description "Compute the displacement of a spring give the load <spring_load> "
  "and spring stiffness <stiffness>";

  local x, out;
  PrintStartProc(procname);

  out := RealDomain[solve](
    spring_load = integrate(spring_stiffness(x), x = 0..Dx), Dx
    );

  PrintEndProc(procname);
  return out;
end proc: # ComputeSpringDisplacement

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeSpringEnergy := proc(
  spring_load::algebraic,      # Load on the spring
  spring_stiffness::procedure, # Spring stiffness
  $)::algebraic;

  description "Compute the potential energy of a spring give the load <spring_load> "
  "and spring stiffness <stiffness>";

  local disp, x, out;
  PrintStartProc(procname);

  disp := ComputeSpringDisplacement(spring_load, spring_stiffness);
  out  := integrate(integrate(spring_stiffness(x), x), x = 0..disp);

  PrintEndProc(procname);
  return out;
end proc: # ComputeSpringEnergy

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeSupportDisplacements := proc(
  obj::SUPPORT, # Support object
  $)::nothing;

  description "Compute the displacements of the support <obj> from its "
    "support reactions";

  local sup_disp, disp_vec, i, disp, x, sup_reac;
  PrintStartProc(procname);

  sup_disp := [];
  disp_vec := [tx, ty, tz, rx, ry, rz];
  sup_reac := [FX, FY, FZ, MX, MY, MZ];

  for i from 1 to 6 do
    if (obj[parse("constrained_dof")][i] = 1) and
        has(obj[parse("support_reactions")], sup_reac[i]) then
      disp := ComputeSpringDisplacement(subs(obj[parse("support_reactions")], - sup_reac[i]),
        (x -> obj[parse("stiffness")](x)[i]));
      sup_disp := sup_disp union [disp_vec[i] = disp];
    end if;
  end do;

  if (verbose_mode > 0) then
    printf(
      "%*sMessage (in ComputeSupportDisplacements) updating %s %s's displacements...\n",
      print_indent, "", obj[parse("type")], obj[parse("name")]
      );
  end if;

  obj[parse("displacements")] := sup_disp;

  if (verbose_mode > 0) then
    printf("%*sDONE\n", print_indent, "");
  end if;

  PrintEndProc(procname);
end proc: # ComputeSupportDisplacements

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
  }, $)::STRUCTURE;

  description "Create a STRUCTURE object with inputs: structure objects <objs>, "
    "external actions <exts>, optional hyperstatic variables <hyper_vars>, optional "
    "hyperstatic displacements <hyper_disp>";

  local num_dof, i, names, candidate_hyp_vars, Graph, out;
  PrintStartProc(procname);

  # Check for duplicate names
  names := [];
  for i from 1 to nops(objs) do
    if member(objs[i][parse("name")], names) then
      error "duplicate names found on structure objects";
    end if;
    names := names union [objs[i][parse("name")]];
  end do;

  num_dof, Graph := ComputeDOF(objs);

  if (num_dof < 0) then
    if (nops(hyper_vars) <> -num_dof) then
      candidate_hyp_vars := [];
      for i from 1 to nops(objs) do
        if IsSupport(objs[i]) or IsJoint(objs[i]) then
          candidate_hyp_vars := [
            op(candidate_hyp_vars),
            op(objs[i][parse("variables")])
            ];
        end if;
      end do;
    WARNING(
      "the structure is hyperstatic with %1 overconstrained directions, "
      "please check the structure supports and joints. Also consider defining "
      "the hyperstatic variables by adding 'hyper_vars' property in the "
      "'MakeStructure' method or simply defining 'hyperstatic_variables' "
      "field in an already existing STRUCTURE object by choosing from the "
      "folloving hyperstatic candidate variables: %2",
      abs(num_dof), candidate_hyp_vars
      );
    else
      if (verbose_mode > 0) then
        printf("%*sMessage (in MakeStructure) "
          "hyperstatic structure detected with %d overconstrained directions\n",
          print_indent, "", abs(num_dof));
      end if;
    end if;
  elif (num_dof > 0) then
    #error "not enough constraints in the structure";
    WARNING(
      "the structure is underconstrained with %1 unconstrained directions. "
      "Results computation may fail due to rigid body motions.",
      num_dof
    );
  else
    if (verbose_mode > 0) then
      printf("%*sMessage (in MakeStructure) isostatic structure detected\n", print_indent, "");
    end if;
  end if;

  out := table({
    parse("type")                      = STRUCTURE,
    parse("objects")                   = objs,
    parse("external_actions")          = exts,
    parse("dof")                       = num_dof,
    parse("connections_graph")         = Graph,
    parse("hyperstatic_variables")     = hyper_vars,
    parse("hyperstatic_displacements") = hyper_disp,
    parse("support_reactions_solved")  = false,
    parse("internal_actions_solved")   = false,
    parse("displacement_solved")       = false
    });

  PrintEndProc(procname);
  return out;
end proc: # MakeStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsStructure := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a STRUCTURE object";

  local out;
  PrintStartProc(procname);

  if (obj[parse("type")] = STRUCTURE) then
    out := true;
  else
    out := false;
  end if;

  PrintEndProc(procname);
  return op(out);
end proc: # IsStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanStructure := proc(
  obj::STRUCTURE, # Object to be cleaned
  $)::nothing;

  description "Clean STRUCTURE object <obj> internal variables";

  local i;
  PrintStartProc(procname);

  # Clean internal variables
  obj[parse("support_reactions_solved")] := false;
  obj[parse("internal_actions_solved")]  := false;
  obj[parse("displacement_solved")]      := false;

  # Clean objects
  for i from 1 to nops(obj[parse("objects")]) do
    if IsBeam(obj[i]) then
      obj[parse("objects")][i] := CleanBeam(i);
    elif IsRod(obj[i]) then
      obj[parse("objects")][i] := CleanRod(i);
    elif IsSupport(obj[i]) then
      obj[parse("objects")][i] := CleanSupport(i);
    elif IsJoint(obj[i]) then
      obj[parse("objects")][i] := CleanJoint(i);
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

  description "Compute the degree of freedom of the input structure objects <objs>";

  local dof, objs_tmp, i, j, k, vertex, colors, G;
  PrintStartProc(procname);

  dof      := 0;
  objs_tmp := objs union [earth];

  # Built connections graph
  vertex := [];
  colors := [];
  if (verbose_mode > 0) then
    printf("%*sMessage (in ComputeDOF) checking structure connections...\n", print_indent, "");
  end if;
  for i from 1 to nops(objs_tmp) do
    vertex := vertex union [objs_tmp[i][parse("name")]];
    colors := colors union [ObjectColor(objs_tmp[i])];
  end do;
  G := GraphTheory[Graph](vertex);
  GraphTheory[HighlightVertex](G, vertex, colors);
  for i from 1 to nops(objs_tmp) do
    if IsSupport(objs_tmp[i]) or IsJoint(objs_tmp[i]) then
      for j from 1 to nops(objs_tmp) do
        if (member(objs_tmp[j][parse("name")], objs_tmp[i][parse("targets")])) then
          GraphTheory[AddEdge](G, {objs_tmp[i][parse("name")], objs_tmp[j][parse("name")]});
        end if;
      end do;
    end if;
  end do;

  # Check graph connections
  if GraphTheory[IsConnected](G) then
    if (verbose_mode > 0) then
      printf("%*sDONE\n", print_indent, "");
    end if;
  else
    error "unconnected elements detected in the structure";
  end if;

  if (verbose_mode > 0) then
    printf("%*sMessage (in ComputeDOF) computing degrees of freedom...\n", print_indent, "");
  end if;

  for i from 1 to nops(objs_tmp) do
    if IsBeam(objs_tmp[i]) then
      dof := dof + 6;
    elif IsRigidBody(objs_tmp[i]) then
      dof := dof + 6;
    elif IsRod(objs_tmp[i]) then
      dof := dof + 5;
    elif IsJoint(objs_tmp[i]) then
      dof := dof - add(objs_tmp[i][constrained_dof][k], k = 1..6) * (nops(objs_tmp[i][parse("targets")]) - 1);
    elif IsSupport(objs_tmp[i]) then
      dof := dof - add(objs_tmp[i][constrained_dof][k], k = 1..6) * (nops(objs_tmp[i][parse("targets")]) - 1);
    end if;
  end do;

  if (verbose_mode > 0) then
    printf("%*sDONE\n", print_indent, "");
    printf("%*sMessage (in ComputeDOF) display degrees of freedom... DOF = %d\n",print_indent, "", dof);
    printf("%*sMessage (in ComputeDOF) display connections graph...\n", print_indent, "");
    if (verbose_mode > 1) then
      print(GraphTheory[DrawGraph](G));
    end if;
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
  obj::STRUCTURE, # Object to be cleaned
  $)::procedure;

  description "Draw the connections graph of the STRUCTURE object <obj>";

  local  out;
  PrintStartProc(procname);

  out := plots:-display(
    GraphTheory[DrawGraph](obj[parse("connections_graph")], layout = tree),
    title = "Structure connections graph");

  PrintEndProc(procname);
  return out;
end proc: # DrawStructureGraph

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

NewtonEuler := proc(
  exts::{ # External actions
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  },
  obj::{BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}, # Object to compute the equilibrium
  {
    pole::list := [0,0,0],                       # Pole to compute the equilibrium
    upper_lim::algebraic := obj[parse("length")] # Upper limit of the integration
  }, $)

  description "Compute the Newton-Euler static equilibrium equations given a set "
    "of external actions <exts>, and object to compute the equilibrium <obj>, "
    "the axial coordinate of the pole <pole>, and an optional upper limit of the "
    "integration <upper_lim>";

  local eq_T, eq_R, i, x, arm, out;
  PrintStartProc(procname);

  eq_T := [0, 0, 0];
  for i from 1 to nops(exts) do
    if exts[i][parse("target")] = obj[parse("name")] then
      if IsForce(exts[i]) then
        eq_T := eq_T + exts[i][parse("components")];
      elif IsQForce(exts[i]) then
        eq_T := eq_T + map(
          integrate, exts[i][parse("components")](x), x = 0..upper_lim
          );
      end if;
    else
      WARNING("Message (in NewtonEuler) %1 is not applied to %2", exts[i], obj);
    end if;
  end do;

  eq_R := [0, 0, 0];
  for i from 1 to nops(exts) do
    if (exts[i][parse("target")] = obj[parse("name")]) then
      if IsMoment(exts[i]) then
        eq_R := eq_R + exts[i][parse("components")];
      elif IsForce(exts[i]) then
          arm := <op(pole - exts[i][parse("coordinate")])>;
        eq_R := eq_R +
          convert(LinearAlgebra[CrossProduct](<op(exts[i][parse("components")])>, arm), list);
      elif IsQForce(exts[i]) then
        arm := <op(pole - [x, 0, 0])>;
        eq_R := eq_R + map(integrate,
          convert(LinearAlgebra[CrossProduct](<op(exts[i][parse("components")](x))>, arm), list),
          x = 0..upper_lim);
      elif IsQMoment(FMQ[i]) then
        eq_R := eq_R + map(integrate,
          exts[i][parse("components")](x), x = 0..upper_lim);
      end if;
    else
      WARNING("Message (in NewtonEuler) %1 is not applied to %2", exts[i], obj);
    end if;
  end do;

  out := eq_T union eq_R;

  PrintEndProc(procname);
  return out;
end proc: # NewtonEuler

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SolveStructure := proc(
  struct::STRUCTURE, # Structure object
  {
    compute_intact::boolean  := false, # Internal actions computation flag
    compute_disp::boolean    := false, # Displacement computation flag
    timoshenko_beam::boolean := false  # Timoshenko beam flag
  }, $)

  description "Solve the static equilibrium of a structure with inputs: "
    "structure <struct>, optional compute internal action enabling flag "
    "<compute_intact>, optional compute displacement enabling flag "
    "<compute_disp>, optional Timoshenko beam flag <timoshenko_beam>";

  local g_load, S_obj, S_rigid, S_ext, S_support, S_joint, S_con_forces, vars, sol, obj,
    x;

  PrintStartProc(procname);

  # Clean structure
  CleanStructure(struct);

  # Parsing inputs
  S_obj        := {};
  S_rigid      := {};
  S_ext        := {};
  S_support    := {};
  S_joint      := {};
  S_con_forces := {};
  vars         := [];
  for obj in struct[parse("objects")] do
    if IsBeam(obj) or IsRod(obj) then
      S_obj := S_obj union {obj};
    elif IsRigidBody(obj) then
      S_rigid := S_rigid union {obj};
    elif IsSupport(obj) then
      S_support    := S_support union {obj};
      S_con_forces := S_con_forces union obj[parse("forces")] union obj[parse("moments")];
      vars         := vars union obj[parse("variables")];
    elif IsJoint(obj) then
      S_joint      := S_joint union {obj};
      S_con_forces := S_con_forces union obj[parse("forces")] union obj[parse("moments")];
      vars         := vars union obj[parse("variables")];
    end if;
      unassign('obj');
  end do;

  S_ext := struct[parse("external_actions")];

  # Add gravity distributed load
  if (gravity <> [0, 0, 0]) then
    for obj in S_obj union S_rigid do
      if IsRod(obj) then
        WARNING("Message (in SolveStructure) gravity load is not supported for rod %1", obj);
      elif IsBeam(obj) then
        g_load||(obj[parse("name")]) := MakeQForce(
          (x -> gravity *~ obj[parse("area")](x) *~ obj[parse("material")][parse("density")]),
          obj,ground
          );
        S_ext := S_ext union {g_load||(obj[parse("name")])};
      elif IsRigidBody(obj) then
        g_load||(obj[parse("name")]) := MakeForce(
          gravity *~ obj[parse("mass")], obj[parse("COM")], obj, ground
          );
        S_ext := S_ext union {g_load||(obj[parse("name")])};
      end if;
    end do;
  end if;

  # Solve isostatic structure
  if (struct[dof] >= 0) then

    if struct[dof] > 0 then
      WARNING("Message (in SolveStructure) structure is underconstrained. "
        "Trying to solve it anyway. Results computation may fail due to rigid "
        "body motions.");
    end if;

    if (verbose_mode > 0) then
      printf("%*sMessage (in SolveStructure) solving the isostatic structure...\n", print_indent, "");
    end if;
    sol := IsostaticSolver(
      S_obj union S_rigid union S_joint union S_support,
      S_ext union S_con_forces,
      vars
      );
    if (verbose_mode > 0) then
      printf("%*sDONE\n", print_indent, "");
      if (verbose_mode > 1) then
        printf("Message (in SolveStructure) solutions:\n");
        print(<sol>);
        printf("%*sMessage (in SolveStructure) updating support reactions fields...\n", print_indent, "");
      end if;
    end if;
    # Update support reactions properties
    for obj in S_support do
      obj[parse("support_reactions")] := [
        seq(lhs(obj[parse("support_reactions")][j]) = subs(sol, rhs(obj[parse("support_reactions")][j])),
        j = 1..nops(obj[parse("support_reactions")]))
      ];
    end do;
    if (verbose_mode > 0) then
      printf("%*sDONE\n", print_indent, "");
    end if;

  # Solve hyperstatic structure
  elif (struct[dof] < 0) then

    if (nops(struct[parse("hyperstatic_variables")]) <> -struct[dof]) then
      error "mismatch in the structure degrees of freedom, check the hyper"
        "static variables of the structure and update the structure object";
    end if;
    if (verbose_mode > 0) then
      printf("%*sMessage (in SolveStructure) solving the hyperstatic structure\n", print_indent, "");
    end if;
    sol := HyperstaticSolver(
      S_obj union S_joint union S_support,
      S_ext union S_con_forces,
      vars,
      struct[parse("hyperstatic_variables")],
      struct[parse("hyperstatic_displacements")],
      parse("timoshenko_beam") = timoshenko_beam
      );
    if (verbose_mode > 0) then
      printf("%*sDONE\n", print_indent, "");
      if (verbose_mode > 1) then
        printf("%*sMessage (in SolveStructure) hyperstatic solver solution:\n", print_indent, "");
        print(<sol>);
      end if;
    end if;
    # Update objects internal actions
    for obj in S_obj do
      obj[parse("internal_actions")] := subs(sol, obj[parse("internal_actions")]);
    end do;
    # Set internal actions computed flag
    struct[parse("internal_actions_solved")] := true;
    # Update support reactions properties
    if (verbose_mode > 0) then
      printf("%*sMessage (in SolveStructure) updating support reactions fields...\n", print_indent, "");
    end if;
    for obj in S_support do
    obj[parse("support_reactions")] := [
      seq(lhs(obj[parse("support_reactions")][j]) = subs(sol,rhs(obj[parse("support_reactions")][j])),
      j = 1..nops(obj[parse("support_reactions")]))
      ];
    end do;
    if (verbose_mode > 0) then
      printf("%*sDONE\n", print_indent, "");
    end if;
  end if;

  # Set support reactions solved flag
  struct[parse("support_reactions_solved")] := true;

  # Compute internal actions
  if ((compute_intact) or (compute_disp)) and not struct[parse("internal_actions_solved")] then
    ComputeInternalActions(
      S_obj, S_ext union S_con_forces, sol
      );
    # Set internal actions computed flag
    struct[parse("internal_actions_solved")] := true;
  end if;

  # Compute displacements
  if (compute_disp) and not struct[parse("displacement_solved")] then
    ComputeDisplacements(
      S_obj union S_joint union S_support, S_ext union S_con_forces, sol
      );
    # Set displacements computed flag
    struct[parse("displacement_solved")] := true;
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
  vars::list,       # Variables
  hyper_vars::list, # Hyperstatic variables
  hyper_disp::list, # Hyperstatic displacements
  {
    timoshenko_beam::boolean := false # Timoshenko beam flag
  }, $)

  description "Solve hyperstatic structure with inputs objects <objs>, external "
    "actions <exts>, variables <vars>, hyperstatic variables <hyper_vars>, "
    "hyperstatic displacements <hyper_disp> and optional Timoshenko beam flag "
    "<timoshenko_beam>";

  local hyper_eq, hyper_load, hyper_comps, hyper_compliant_disp, hyper_support, i, obj,
        iso_vars, iso_sol, hyper_sol, sol, P, S_objs, E_objs, out;
  PrintStartProc(procname);

  # Parse input objects and find objects with internal actions property
  S_objs := [seq(
    `if`(IsBeam(objs[i]) or IsRod(objs[i]), objs[i], NULL),
    i = 1..nops(objs))
    ];

  if (verbose_mode > 1) then
    printf("%*sMessage (in HyperstaticSolver) solving the hyperstatic variables...\n", print_indent, "");
  end if;

  # Create a solution as function of the hyperstatic variables
  iso_vars := [seq(
    `if`(member(vars[i], hyper_vars), NULL, vars[i]),
    i = 1..nops(vars))
    ];
  iso_sol := IsostaticSolver(objs, exts, iso_vars);

  # Compute internal actions
  ComputeInternalActions(S_objs, exts, iso_sol);

  # Extract the deformable objects
  E_objs := [seq(
    `if`(not IsRigidBody(objs[i]), objs[i], NULL),
    i = 1..nops(objs))
    ];

  # Compute structure internal energy
  P := ComputePotentialEnergy(E_objs, iso_sol, parse("timoshenko_beam") = timoshenko_beam);

  hyper_eq := [];

  for i from 1 to nops(hyper_vars) do
    # Compose the hyperstatic equation
    hyper_eq := hyper_eq union [diff(P, hyper_vars[i]) = hyper_disp[i]];
  end do;

  # Solve hyperstatic equations
  hyper_sol := op(RealDomain[solve](convert(hyper_eq, signum), hyper_vars));

  if (verbose_mode > 1) then
    printf("%*sDONE\n", print_indent, "");
  end if;

  out := hyper_sol union subs(hyper_sol, iso_sol);

  PrintEndProc(procname);
  return out;
end proc: # HyperstaticSolver

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputePotentialEnergy := proc(
  objs::{ # Structure objects
    list({BEAM, ROD, SUPPORT, JOINT}),
    set( {BEAM, ROD, SUPPORT, JOINT})
  },
  sol::{list(`=`),set(`=`)} := [], # Substitutions
  {
    timoshenko_beam := false # Timoshenko beam flag
  }, $)

  description "Compute the internal potential energy of the structure given the "
    "objects <objs> and optional Timoshenko beam flag <timoshenko_beam>";

  local obj, P, x;
  PrintStartProc(procname);

  P := 0;
  for obj in objs do
    if IsBeam(obj) or IsRod(obj) then
      # Normal action N contribution
      if (member(N, map(lhs, obj[parse("internal_actions")]))) and
          (subs(obj[parse("internal_actions")](x), N(x)) <> 0) then
        P := P + integrate(
          subs(obj[parse("internal_actions")](x),
            N(x)^2/(2*obj[parse("material")][parse("elastic_modulus")]*obj[parse("area")](x))
          ), x = 0..obj[parse("length")]);
      end if;
      if timoshenko_beam then
        # Shear action Ty contribution
        if (member(Ty, map(lhs, obj[parse("internal_actions")]))) and
            (subs(obj[parse("internal_actions")](x), Ty(x)) <> 0) then
          P := P + integrate(
            subs(obj[parse("internal_actions")](x),
              Ty(x)^2/(2*obj[parse("timo_shear_coeff")](x)[1]*obj[parse("material")][parse("shear_modulus")]*obj[parse("area")](x))
            ), x = 0..obj[parse("length")]);
        end if;
        # Shear action Tz contribution
        if (member(Tz, map(lhs, obj[parse("internal_actions")]))) and
            (subs(obj[parse("internal_actions")](x), Tz(x)) <> 0) then
          P := P + integrate(
            subs(obj[parse("internal_actions")](x),
              Tz(x)^2/(2*obj[parse("timo_shear_coeff")](x)[2]*obj[parse("material")][parse("shear_modulus")]*obj[parse("area")](x))
            ), x = 0..obj[parse("length")]);
        end if;
      end if;
      # Bending moment action Mx contribution
      if (member(Mx, map(lhs, obj[parse("internal_actions")]))) and
          (subs(obj[parse("internal_actions")](x), Mx(x)) <> 0) then
        P := P + integrate(
          subs(obj[parse("internal_actions")](x),
            Mx(x)^2/(2*obj[parse("material")][parse("shear_modulus")]*obj[parse("inertias")][1](x))
          ), x = 0..obj[parse("length")]);
          end if;
      # Bending moment action My contribution
      if (member(My, map(lhs, obj[parse("internal_actions")]))) and
          (subs(obj[parse("internal_actions")](x), My(x)) <> 0) then
        P := P + integrate(
          subs(obj[parse("internal_actions")](x),
            My(x)^2/(2*obj[parse("material")][parse("elastic_modulus")]*obj[parse("inertias")][3](x))
          ), x = 0..obj[parse("length")]);
      end if;
      # Bending moment action Mz contribution
      if (member(Mz, map(lhs, obj[parse("internal_actions")]))) and
          (subs(obj[parse("internal_actions")](x), Mz(x)) <> 0) then
        P := P + integrate(
          subs(obj[parse("internal_actions")](x),
            Mz(x)^2/(2*obj[parse("material")][parse("elastic_modulus")]*obj[parse("inertias")][2](x))
          ), x = 0..obj[parse("length")]);
      end if;
    elif IsCompliantSupport(obj) then
      # Support reaction Fx contribution
      if (subs(obj[parse("support_reactions")], FX) <> 0) and
          (obj[parse("stiffness")](x)[1] <> infinity) and
           (obj[parse("constrained_dof")][1] <> 0) then
        P := P + subs(obj[parse("support_reactions")], sol, ComputeSpringEnergy(-FX, (x -> obj[parse("stiffness")](x)[1])));
      end if;
      # Support reaction Fy contribution
      if (subs(obj[parse("support_reactions")], FY) <> 0) and
          (obj[parse("stiffness")](x)[2] <> infinity) and
           (obj[parse("constrained_dof")][2] <> 0) then
        P := P + subs(obj[parse("support_reactions")], sol, ComputeSpringEnergy(-FY, (x -> obj[parse("stiffness")](x)[2])));
      end if;
      # Support reaction Fz contribution
      if (subs(obj[parse("support_reactions")], FZ) <> 0) and
          (obj[parse("stiffness")](x)[3] <> infinity) and
           (obj[parse("constrained_dof")][3] <> 0) then
        P := P + subs(obj[parse("support_reactions")], sol, ComputeSpringEnergy(-FZ, (x -> obj[parse("stiffness")](x)[3])));
      end if;
      # Support reaction Mx contribution
      if (subs(obj[parse("support_reactions")], MX) <> 0) and
          (obj[parse("stiffness")](x)[4] <> infinity) and
           (obj[parse("constrained_dof")][4] <> 0) then
        P := P + subs(obj[parse("support_reactions")], sol, ComputeSpringEnergy(-MX, (x -> obj[parse("stiffness")](x)[4])));
      end if;
      # Support reaction My contribution
      if (subs(obj[parse("support_reactions")], MY) <> 0) and
          (obj[parse("stiffness")](x)[5] <> infinity) and
           (obj[parse("constrained_dof")][5] <> 0) then
        P := P + subs(obj[parse("support_reactions")], sol, ComputeSpringEnergy(-MY, (x -> obj[parse("stiffness")](x)[5])));
      end if;
      # Support reaction Mz contribution
      if (subs(obj[parse("support_reactions")], MZ) <> 0) and
          (obj[parse("stiffness")](x)[6] <> infinity) and
           (obj[parse("constrained_dof")][6] <> 0) then
        P := P + subs(obj[parse("support_reactions")], sol, ComputeSpringEnergy(-MZ, (x -> obj[parse("stiffness")](x)[6])));
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
  vars::list, # Variables to solve
  $)

  description "Solve the isostatic structure equilibrium equation system given "
    "the structure objects <objs>, the external actions <exts> and the variables "
    "<vars> to solve";

  local eq, i, j, active_ext, sol, A, B, rank_eq, vars_tmp;
  PrintStartProc(procname);

  # Compute structure equations
  if (verbose_mode > 1) then
    printf("%*sMessage (in IsostaticSolver) computing the equilibrium equation for "
      "the isostatic structure...\n", print_indent, "");
  end if;
  eq := [];
  for i from 1 to nops(objs) do
    active_ext := {};
    for j from 1 to nops(exts) do
      if (exts[j][parse("target")] = objs[i][parse("name")]) then
        active_ext := active_ext union {exts[j]};
      end if;
    end do;
    eq := eq union NewtonEuler(active_ext, objs[i]);

    # Add joints and supports constraint equations
    if IsSupport(objs[i]) or IsJoint(objs[i]) then
      eq := eq union objs[i][parse("constraint_loads")];
    end if;
  end do;

  # Remove NULL equations
  eq := remove(x -> x = 0, simplify(eq));

  # Remove non used variables
  vars_tmp := vars;
  for i from 1 to nops(vars) do
  if (has(eq, vars[i])) = false then
    vars_tmp := remove(x -> x = vars[i], vars_tmp);
    WARNING(
      "Message (in IsostaticSolver) %1 was removed from variables because it "
      "is not used in the equations",
      vars[i]
      );
      end if;
  end do;
  if (verbose_mode > 1) then
    printf("%*sDONE\n", print_indent, "");
  end if;

  # Structure equations check
  A, B    := LinearAlgebra[GenerateMatrix](eq, vars_tmp);
  rank_eq := LinearAlgebra[Rank](A);

  if (verbose_mode > 1) then
    printf("%*sMessage (in IsostaticSolver) structure equilibrium equations:\n", print_indent, "");
    print(<op(eq)>);
    printf("%*sMessage (in IsostaticSolver) structure unknown variables:\n", print_indent, "");
    print(vars_tmp);
  end if;

  if (rank_eq <> nops(vars_tmp)) or (nops(eq) <> nops(vars_tmp)) then
    error "inconsistent system of equation, got %1 equations and %2 variables. "
      "Rank of the system  wrt the system variables is %3. Check structure "
      "supports and joints",
      nops(eq), nops(vars_tmp), rank_eq;
  end if;

  # Solve structure equations
  if (verbose_mode > 1) then
    printf("%*sMessage (in IsostaticSolver) computing the structure reaction forces...\n", print_indent, "");
  end if;

  sol := simplify(op(RealDomain[solve](eq, vars_tmp)));

  if (verbose_mode > 1) then
    printf("%*sDONE\n", print_indent, "");
  end if;

  PrintEndProc(procname);
  return sol;
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
  $)::nothing;

  description "Programmatic computation of internal actions for structure"
    "objects with given external actions and structure solution";

  local i, j, active_ext, subs_ext;

  PrintStartProc(procname);

  # Substitute structure solution into loads
  subs_ext := map2(subs, sol, map(op,exts));

  for i from 1 to nops(objs) do
    # Extract active loads
    active_ext := {};
    for j from 1 to nops(subs_ext) do
      if (subs_ext[j][parse("target")] = objs[i][parse("name")]) then
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
  }, $)::nothing;

  description "Programmatic computation of internal actions for structure "
    "object with given external actions, it returns the internal actions as "
    "function of the axial variable 'x'";

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
  # computation, but results must be considered valid only in the assumed range
   Physics[Assume](x > 0, x < obj[parse("length")]);

  # Compute internal actions for concentrated loads as effect overlay
  for i from 1 to nops(exts) do
    if IsForce(exts[i]) then
      N_sol  := `simplify/piecewise`(N_sol  - piecewise(x >= exts[i][parse("coordinate")][1] and x <= obj[parse("length")], exts[i][parse("components")][1]), x);
      Ty_sol := `simplify/piecewise`(Ty_sol + piecewise(x >= exts[i][parse("coordinate")][1] and x <= obj[parse("length")], exts[i][parse("components")][2]), x);
      Tz_sol := `simplify/piecewise`(Tz_sol + piecewise(x >= exts[i][parse("coordinate")][1] and x <= obj[parse("length")], exts[i][parse("components")][3]), x);
      My_sol := `simplify/piecewise`(My_sol + integrate(piecewise(x >= exts[i][parse("coordinate")][1] and x <= obj[parse("length")], exts[i][parse("components")][3]), x = 0..x), x);
      Mz_sol := `simplify/piecewise`(Mz_sol + integrate(piecewise(x >= exts[i][parse("coordinate")][1] and x <= obj[parse("length")], exts[i][parse("components")][2]), x = 0..x), x);
    elif IsMoment(exts[i]) then
      Mx_sol := `simplify/piecewise`(Mx_sol - piecewise(x >= exts[i][parse("coordinate")][1] and x <= obj[parse("length")], exts[i][parse("components")][1]), x);
      My_sol := `simplify/piecewise`(My_sol + piecewise(x >= exts[i][parse("coordinate")][1] and x <= obj[parse("length")], exts[i][parse("components")][2]), x);
      Mz_sol := `simplify/piecewise`(Mz_sol - piecewise(x >= exts[i][parse("coordinate")][1] and x <= obj[parse("length")], exts[i][parse("components")][3]), x);
    elif IsQForce(exts[i]) then
      N_sol  := `simplify/piecewise`(N_sol  - integrate(exts[i][parse("components")](x)[1], x = 0..x), x);
      Ty_sol := `simplify/piecewise`(Ty_sol + integrate(exts[i][parse("components")](x)[2], x = 0..x), x);
      Tz_sol := `simplify/piecewise`(Tz_sol + integrate(exts[i][parse("components")](x)[3], x = 0..x), x);
      My_sol := `simplify/piecewise`(My_sol + integrate(integrate(exts[i][parse("components")](x)[3], x = 0..x), x = 0..x), x);
      Mz_sol := `simplify/piecewise`(Mz_sol + integrate(integrate(exts[i][parse("components")](x)[2], x = 0..x), x = 0..x), x);
    elif IsQMoment(FMQ[i]) then
      Mx_sol := `simplify/piecewise`(Mx_sol - integrate(exts[i][parse("components")](x)[1], x = 0..x), x);
      My_sol := `simplify/piecewise`(My_sol + integrate(exts[i][parse("components")](x)[2], x = 0..x), x);
      Mz_sol := `simplify/piecewise`(Mz_sol - integrate(exts[i][parse("components")](x)[3], x = 0..x), x);
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
    print_indent, "", obj[parse("type")], obj[parse("name")]
    );
  end if;

  obj[parse("internal_actions")] := ia;

  if (verbose_mode > 1) then
    printf("%*sDONE\n", print_indent, "");
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
    timoshenko_beam::boolean := false # Timoshenko beam flag
  }, $)::nothing;

  description "Compute the Structure displacements";

  local obj, x, X, disp, rx_sol, ry_sol, rz_sol, ux_sol, uy_sol, uz_sol;
  PrintStartProc(procname);

  # Cicle on the structure objects
  for obj in objs do
    # Beam
    if IsBeam(obj) then
      # Compute displacements
      rx_sol := integrate(subs(obj[parse("internal_actions")](x), Mx(x)/(obj[parse("material")][parse("shear_modulus")]*obj[parse("inertias")][1](x))), x = 0..X);
      ry_sol := integrate(subs(obj[parse("internal_actions")](x), My(x)/(obj[parse("material")][parse("elastic_modulus")]*obj[parse("inertias")][2](x))), x = 0..X);
      rz_sol := integrate(subs(obj[parse("internal_actions")](x), Mz(x)/(obj[parse("material")][parse("elastic_modulus")]*obj[parse("inertias")][3](x))), x = 0..X);
      ux_sol := integrate(subs(obj[parse("internal_actions")](x), N(x)/(obj[parse("material")][parse("elastic_modulus")]*obj[parse("area")](x))), x = 0..X);
      uy_sol := integrate(eval(ry_sol, X = x), x = 0..X);
      uz_sol := integrate(eval(rz_sol, X = x), x = 0..X);
      if timoshenko_beam then
        uy_sol := ty_sol + integrate(subs(obj[parse("internal_actions")](x), Ty(x)/(obj[parse("timo_shear_coeff")](x)[1]*obj[parse("material")][parse("shear_modulus")]*obj[parse("area")](x))), x = 0..X);
        uz_sol := tz_sol + integrate(subs(obj[parse("internal_actions")](x), Tz(x)/(obj[parse("timo_shear_coeff")](x)[2]*obj[parse("material")][parse("shear_modulus")]*obj[parse("area")](x))), x = 0..X);
      end if;
      disp := [
        rx = unapply(rx_sol, X), ry = unapply(ry_sol, X), rz = unapply(rz_sol, X),
        ux = unapply(ux_sol, X), uy = unapply(uy_sol, X), uz = unapply(uz_sol, X)
      ];

      # Update object displacements
      obj[parse("displacements")] := disp;

    # Rod
    elif IsRod(obj) then
      # Compute displacements
      ux_sol := integrate(subs(obj[parse("internal_actions")](x), N(x)/(obj[parse("material")][parse("elastic_modulus")]*obj[parse("area")](x))), x = 0..X);
      disp := [
        ux = unapply(ux_sol, X)
      ];

      # Update object displacements
      obj[parse("displacements")] := disp;

    # Support
    elif IsCompliantSupport(obj) then
      # Compute displacements
      ComputeSupportDisplacements(obj);
    end if;
  end do;

  PrintEndProc(procname);
end proc: # ComputeDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ObjectColor := proc(
  obj::{BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH}, # Object to be colored
$)::string;

description "Return the color of the object <obj>";

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
  elif IsJoint(obj) then
    color := Joint_color;
  elif IsEarth(obj) then
    color := Earth_color;
  end if;

  PrintEndProc(procname);
  return color;
end proc: # ObjectColor


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotBeam := proc(
  obj::BEAM, # Beam to be plot
  {
    data::{list(`=`),set(`=`)} := [] # Substitutions
  },
  $)::procedure;

  description "Plot a the SUPPORT object <obj>";

  local P1, P2, out;
  PrintStartProc(procname);

  P1 := subs(data, Origin(obj[parse("frame")]));
  P2 := subs(data, Origin(obj[parse("frame")].Translate(obj[parse("length")], 0, 0)));

  out := plots:-display(
    plottools:-line(convert(P1[1..3], list), convert(P2[1..3], list), thickness = 6),
    linestyle = solid, color = ObjectColor(obj));

  PrintEndProc(procname);
  return out;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotRod := proc(
  obj::ROD, # Rod to be plot
  {
    data::{list(`=`),set(`=`)} := [] # Substitutions
  },
  $)::procedure;

  description "Plot a the ROD object <obj>";

  local P1, P2, out;
  PrintStartProc(procname);

  P1 := subs(data, Origin(obj[parse("frame")]));
  P2 := subs(data, Origin(obj[parse("frame")].Translate(obj[parse("length")], 0, 0)));

  out := plots:-display(
    plottools:-line(convert(P1[1..3], list), convert(P2[1..3], list), thickness = 4),
    linestyle = dash, color = ObjectColor(obj));

  PrintEndProc(procname);
  return out;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotJoint := proc(
  obj::JOINT,            # Joint to be plot
  {
    data::{list(`=`),set(`=`)} := [] # Substitutions
  },
  $)::procedure;

  description "Plot a the JOINT object <obj>";

  local O, out;
  PrintStartProc(procname);

  O := subs(data, Origin(
    parse(obj[parse("targets")][1])[parse("frame")].
    Translate(op(ListPadding(obj[parse("coordinates")][1],3)))
    ));

  out := plots:-display(
    plottools:-point(convert(O[1..3], list), symbol='solidsphere', symbolsize = 20),
    linestyle = solid, color = ObjectColor(obj));

  PrintEndProc(procname);
  return out;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotSupport := proc(
  obj::SUPPORT, # Joint to be plot
  {
    data::{list(`=`),set(`=`)} := [] # Substitutions
  },
  $)::procedure;

  local O, out;
  PrintStartProc(procname);

  if nops(obj[parse("targets")])>1 then
    O := subs(data, Origin(
      parse(obj[parse("targets")][2])[parse("frame")].
      Translate(op(ListPadding(obj[parse("coordinates")][2],3)))
      ));
  else
    O := subs(data, Origin(
      earth[parse("frame")].
      Translate(op(ListPadding(obj[parse("coordinates")][1],3)))
      ));
  end if;

  out := plots:-display(
    plottools:-point(convert(O[1..3], list), symbol='solidbox', symbolsize = 20),
    linestyle = solid, color = ObjectColor(obj));

  PrintEndProc(procname);
  return out;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotStructure := proc(
  str::STRUCTURE,                   # Structure to be plot
  data::{list(`=`),set(`=`)} := [], # Substitutions
  $)::list(procedure);

  description "Plot a the STRUCTURE object <obj>";

  local out, obj;
  PrintStartProc(procname);

  out := []:
  for obj in str[parse("objects")] do
    if IsBeam(obj) then
      out := [op(out), PlotBeam(obj, parse("data") = data)];
    elif IsRod(obj) then
      out := [op(out), PlotRod(obj, parse("data") = data)];
    elif IsSupport(obj) then
      out := [op(out), PlotSupport(obj, parse("data") = data)];
    elif IsJoint(obj) then
      out := [op(out), PlotJoint(obj, parse("data") = data)];
    end if;
  end do;

  PrintEndProc(procname);
  return out;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsInsideJoint := proc(
  obj::JOINT,        # Joint object
  p::list,           # Point to be checked
  tol::real := 1e-3, # Tolerance
  $)::boolean;

  description "Check if the point <p> is inside the JOINT <obj>";

  local O, out;
  PrintStartProc(procname);

  if not (nops(p) = 3) then
    error "The input point must be a list of 3 elements";
  end if:

  if (nops(obj[parse("targets")]) > 1) then
    O := Origin(
      parse(obj[parse("targets")][2])[parse("frame")].
      Translate(obj[parse("coordinates")][2], 0, 0)
      );
  else
    WARNING("The support has no targets");
  end if;

  out := (norm(p - O) <= tol);

  PrintEndProc(procname);
  return out;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsInsideSupport := proc(
  obj::SUPPORT,      # Support object
  p::list,           # Point to be checked
  tol::real := 1e-3, # Tolerance
  $)::boolean;

  description "Check if the point <p> is inside the SUPPORT <obj>";

  local O, out;
  PrintStartProc(procname);

  if not (nops(p) = 3) then
    error "The input point must be a list of 3 elements";
  end if:

  if (nops(obj[parse("targets")]) > 1) then
    O := Origin(
      parse(obj[parse("targets")][2])[parse("frame")].
      Translate(obj[parse("coordinates")][2], 0, 0)
      );
  else
    WARNING("The support has no targets");
  end if;

  out := (norm(p - O) <= tol);

  PrintEndProc(procname);
  return out;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsInsideRod := proc(
  obj::ROD, # Rod object
  p::list,  # Point to be checked
  $)::boolean;

  description "Check if the point <p> is inside the ROD <obj>";

  local O, V, W, out;
  PrintStartProc(procname);

  if not (nops(p) = 3) then
    error "The input point must be a list of 3 elements";
  end if:

  O := Origin(obj[parse("frame")]);
  V := obj[parse("frame")].Translate(obj[parse("length")], 0, 0) - O;
  W := p - O;

  out := (dot(W, V) >= 0) and (dot(W, V) <= dot(V, V));

  PrintEndProc(procname);
  return out;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsInsideBeam := proc(
  obj::BEAM, # Beam object
  p::list,   # Point to be checked
  $)::boolean;

  description "Check if the point <p> is inside the BEAM <obj>";

  local O, V, W, out;
  PrintStartProc(procname);

  if not (nops(p) = 3) then
    error "The input point must be a list of 3 elements";
  end if:

  O := Origin(obj[parse("frame")]);
  V := obj[parse("frame")].Translate(obj[parse("length")], 0, 0) - O;
  W := p - O;

  out := (dot(W, V) >= 0) and (dot(W, V) <= dot(V, V));

  PrintEndProc(procname);
  return out;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsInsideStructure := proc(
  obj::STRUCTURE, # Structure object
  p::list,        # Point to be checked
  $)::boolean;

  description "Check if the point <p> is inside the STRUCTURE <obj>";

  local i, out;
  PrintStartProc(procname);

  out := false;
  for i in str[parse("objects")] do
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
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module:

# That's all folks!
