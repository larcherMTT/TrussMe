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
        IsEarth,
        Show,
        Rotate,
        Translate,
        Project,
        InverseFrame,
        IsFrame,
        IsVector,
        IsPoint,
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
        PlotStructure;

global  _gravity,
        ground;

local   ModuleLoad,
        ModuleUnload,
        EARTH,
        GetNames,
        CleanJoint,
        CleanSupport,
        CleanRod,
        CleanBeam,
        CleanStructure,
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
        PlotBeam,
        PlotRod,
        PlotJoint,
        PlotSupport,
        lib_base_path,
        print_indent,
        print_increment;

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
  TypeTools[AddType](FRAME, IsFrame);
  TypeTools[AddType](VECTOR, IsVector);
  TypeTools[AddType](POINT, IsPoint);
  TypeTools[AddType](EARTH, IsEarth);
  TypeTools[AddType](BEAM, IsBeam);
  TypeTools[AddType](ROD, IsRod);
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
             <0, 0, 0, 1>>:

  _gravity := [0, 0, 0]:

  print_indent := 0:
  print_increment := 4:

  EARTH := table({
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
    'VECTOR',
    'POINT',
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

#    ___                       _                 
#   / _ \ _ __   ___ _ __ __ _| |_ ___  _ __ ___ 
#  | | | | '_ \ / _ \ '__/ _` | __/ _ \| '__/ __|
#  | |_| | |_) |  __/ | | (_| | || (_) | |  \__ \
#   \___/| .__/ \___|_|  \__,_|\__\___/|_|  |___/
#      |_|                                     

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

`union` := proc(
  A::{list, set}, # Object A to be united
  B::{list, set},  # Object B to be united
  $) option overload;

  description "Extension of union operator to list objects <A> and <B>";

  if type(A, 'set') and type(B, 'set') then
    return {op(A), op(B)};
  else
    return [op(A), op(B)];
  end if:
end proc: # union

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#   _____                 _   _
#  |  ___|   _ _ __   ___| |_(_) ___  _ __  ___
#  | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
#  |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
#  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsEarth := proc(
  obj::anything, # Object to be tested
  $)::boolean;

  description "Test if an object <obj> is the EARTH object";

  if type(obj, table) and
     (obj[parse("type")] = EARTH) and
     (obj[parse("length")] = 0) and
     (obj[parse("frame")] = ground) and
     (obj[parse("admissible_loads")] = [1, 1, 1, 1, 1, 1]) then
    return true;
  else
    return false;
  end if:
end proc: # IsEarth

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Norm2 := proc(
  v::{list, vector}, # Vector for which the norm is computed
  $)::algebraic;

  description "Compute the norm of a vector <v>";

  return sqrt(add(x, x in v^~2));
end proc: # Norm2

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Show := proc(
  tab::table, # Table to be shown
  $)::nothing;

  description "Show the content of a table <tab>";

  print(tab = tab[parse("type")](op(op(tab))));
end proc: # Show

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

GetNames := proc(
  objs::{ # Structural elements
    list({MATERIAL, BEAM, ROD, SUPPORT, JOINT, EARTH}),
    set( {MATERIAL, BEAM, ROD, SUPPORT, JOINT, EARTH})
  }, $)::{list({string}), set({string})};

  description "Get names of a list/set of objects <objs>";

  if type(objs, 'set') then
    return {seq(objs[i][parse("name")], i = 1..nops(objs))};
  else
    return [seq(objs[i][parse("name")], i = 1..nops(objs))];
  end if:
end proc: # GetNames

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

InverseFrame := proc(
  RF::FRAME, # Reference frame to be inverted
  $)::FRAME;

  description "Inverse transformation matrix of <RF>";

  return LinearAlgebra:-Transpose(RF);
end proc: # InverseFrame

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsFrame := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the input object <obj> is a FRAME object";

  if (type(obj, 'Matrix')) and
     (LinearAlgebra:-RowDimension(obj) = 4) and
     (LinearAlgebra:-ColumnDimension(obj) = 4) and
     (simplify(combine(LinearAlgebra:-Determinant(obj) = 1))) then
    return true;
  else
    return false;
  end if;
end proc: # IsFrame

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsVector := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the input object <obj> is a vector";

  if (type(obj, 'Matrix')) and
     (LinearAlgebra:-Dimension(obj) = 4) and
     (obj[4] = 0) then
    return true;
  else
    return false;
  end if:
end proc: # IsVector

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsPoint := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the input object <obj> is a point";

  if (type(obj, 'Matrix')) and
     (LinearAlgebra:-Dimension(obj) = 4) and
     (obj[4] = 1) then
    return true;
  else
    return false;
  end if:
end proc: # IsPoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Rotate := proc(
  axis::symbol,     # Rotation axis
  angle::algebraic, # Rotation angle (rad)
  $)::FRAME;

  description "Transformation matrix corresponding to the rotation <angle> "
    "around the given <axis>";

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
    error "invalid axis detected";
  end
end proc: # Rotate

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Translate := proc(
  x::algebraic, # X-axis translation component
  y::algebraic, # Y-axis translation component
  z::algebraic, # Z-axis translation component
  $)::FRAME;

  description "Transformation matrix corresponding to the translation <x,y,z>";

  return <<1, 0, 0, 0>|
          <0, 1, 0, 0>|
          <0, 0, 1, 0>|
          <x, y, z, 1>>;
end proc: # Translate

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Origin := proc(
  RF::FRAME, # Reference frame
  $)::vector;

  description "Extract the origin of the reference frame <RF>";

  return <RF[1,4], RF[2,4], RF[3,4], 1>;
end proc: # Translate

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Uvec := proc(
  axis::symbol,        # Axis of the unit vector
  RF::FRAME := ground, # Reference frame
  $)::vector;

  description "Extract the unit vector of the reference frame <RF> along the "
    "given <axis>";

  if (axis = 'X') then
    return <RF[1,1], RF[2,1], RF[3,1], 0>;
  elif (axis = 'Y') then
    return <RF[1,2], RF[2,2], RF[3,2], 0>;
  elif (axis = 'Z') then
    return <RF[1,3], RF[2,3], RF[3,3], 0>;
  else
    error "invalid axis detected";
  end if:
end proc: # Uvec

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

UvecX := proc(
  RF::FRAME := ground, # Reference frame
  $)::vector;

  description "Extract the x-axis unit vector of the reference frame <RF>";

  return <RF[1,1], RF[2,1], RF[3,1], 0>;
end proc: # UvecX

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

UvecY := proc(
  RF::FRAME := ground, # Reference frame
  $)::vector;

  description "Extract the y-axis unit vector of the reference frame <RF>";

  return <RF[1,2], RF[2,2], RF[3,2], 0>;
end proc: # UvecY

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

UvecZ := proc(
  RF::FRAME := ground, # Reference frame
  $)::vector;

  description "Extract the z-axis unit vector of the reference frame <RF>";

  return <RF[1,3], RF[2,3], RF[3,3], 0>;
end proc: # UvecZ

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Project := proc(
  x::{list, vector}, # Vector/point to be projected
  RF_ini::FRAME,     # Reference frame from which the vector/point is expressed
  RF_end::FRAME,     # Reference frame to which the vector/point will be expressed
  $)::{list, vector};

  description "Project <x,y,z>, or vector <x,y,z,0>, or point <x,y,z,1> from "
    "reference frame <RF_ini> to reference frame <RF_end>";

  local x_tmp;

  # Pad input vector with 0 if its length is 3
  if (nops(x) = 3) then
    x_tmp := <x[1], x[2], x[3], 0>;
  elif (nops(x) = 4) then
    x_tmp := <x[1], x[2], x[3], x[4]>;
  else
    error "invalid input vector/point <x> detected";
  end if;

  LinearAlgebra[MatrixInverse](RF_end).RF_ini.x_tmp;
  return simplify([seq(%[i], i = 1..nops(x))]);
end proc: # Project

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeMaterial := proc({
    name::string            := "steel",      # Name of the material
    elastic_modulus::algebraic := 210.0E+09, # Elastic modulus (Pa)
    poisson_modulus::algebraic := 0.3,       # Poisson modulus (-)
    shear_modulus::algebraic   := 77.0E+09,  # Shear modulus (Pa)
    density::algebraic         := 7.4E+03    # Density (kg/m^3)
  }, $)::MATERIAL;

  description "Define a MATERIAL object with inputs: name of the material, "
    "elastic modulus (default = 210.0E9 Pa), Poisson modulus (default = "
    "0.3), shear modulus (default = E/(2*(1+nu))), density (default = "
    "7.4E3 kg/m^3)";

  return table({
    parse("type")            = MATERIAL,
    parse("name")            = name,
    parse("elastic_modulus") = elastic_modulus,
    parse("poisson_modulus") = poisson_modulus,
    parse("shear_modulus")   = shear_modulus,
    parse("density")         = density
    });
end proc: # DefineMaterial

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsMaterial := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the input object <obj> is a MATERIAL object";

  if (obj[parse("type")] = MATERIAL) and
     type(obj, table) and
     type(obj[parse("name")], string) and
     type(obj[parse("elastic_modulus")], algebraic) and
     type(obj[parse("poisson_modulus")], algebraic) and
     type(obj[parse("shear_modulus")], algebraic) and
     type(obj[parse("density")], algebraic) then
    return true;
  else
    return false;
  end if;
end proc: # IsMaterial

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeForce := proc(
  components::list,                        # Force components
  ell::algebraic,                          # Axial coordinate
  obj::{BEAM, ROD, SUPPORT, JOINT, EARTH}, # Target object
  RF::FRAME := ground,                     # Reference frame
  $)::FORCE;

description "Define a FORCE object with inputs: force components, force "
    "application axial coordinate <ell> in [0,L], target object <obj>, optional "
    "reference frame <RF> in which the force is defined (default = ground)";

  local proj_components;

  if IsBeam(obj) or IsRod(obj) then
    if (not type(indets(ell),set(symbol))) then
      if (evalf(ell) < 0) or (evalf(ell) > evalf(obj[parse("length")])) then
        error "force application point must be in [0,L] range";
      end if;
    end if;
  end if;

  proj_components := Project(components, RF, obj[parse("frame")]);
  if IsRod(obj) then
    if (proj_components[2] <> 0) or (proj_components[3] <> 0) then
      error "only axial forces are accepted in ROD objects";
    end if;
  elif IsSupport(obj) or IsJoint(obj) then
    if (evalf(ell) <> 0) then
      error "only null axial coordinate is accepted for SUPPORT and JOINT "
        "objects";
    end if;
  end if;

  return table({
    parse("type")       = FORCE,
    parse("components") = proj_components,
    parse("coordinate") = ell,
    parse("target")     = obj[parse("name")]
    });
end proc: # MakeForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsForce := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a FORCE object";

  if (obj[parse("type")] = FORCE) and
     type(obj, table) and
     type(obj[parse("components")], list) and
     type(obj[parse("coordinate")], algebraic) and
     type(obj[parse("target")], string) then
    return true;
  else
    return false;
  end if;
end proc: # IsForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeMoment := proc(
  components::list,                   # Moment components
  ell::algebraic,                     # Axial coordinate
  obj::{BEAM, SUPPORT, JOINT, EARTH}, # Target object
  RF::FRAME := ground,                # Reference frame
  $)::MOMENT;

  description "Define a MOMENT object with inputs: moment components, "
    "moment application axial coordinate [0,L], target object, optional "
    "reference frame in which the moment is defined (default = ground)";

  local proj_components;

  # FIXME: consider the case of symbolic length or ell
  if IsBeam(obj) then
    if (not type(indets(ell),set(symbol))) then
      if (evalf(ell) < 0) or (evalf(ell) > evalf(obj[parse("length")])) then
        error "moment application point must be in [0,L] range";
      end if;
    end if;
  end if;

  proj_components := Project(components, RF, obj[parse("frame")]);
  if IsSupport(obj) or IsJoint(obj) then
    if (evalf(ell) <> 0) then
      error "only null axial coordinate is accepted for SUPPORT and JOINT "
        "objects";
    end if;
  end if;

  return table({
    parse("type")       = MOMENT,
    parse("components") = proj_components,
    parse("coordinate") = ell,
    parse("target")     = obj[parse("name")]
    });
end proc: # MakeMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsMoment := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a MOMENT object";

  if (obj[parse("type")] = MOMENT) and
     type(obj, table) and
     type(obj[parse("components")], list) and
     type(obj[parse("coordinate")], algebraic) and
     type(obj[parse("target")], string) then
    return true;
  else
    return false;
  end if;
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

  description "Define a 'QFORCE' object with inputs: distributed load components, "
    "target object, initial and final application points (axial coordinates), "
    "optional reference frame in which the load components are defined "
    "(default = ground)";

  local proj_components, x;

  if type(components, procedure) then
    proj_components := unapply(Project(components(x), RF, obj[parse("frame")]),x);
  else
    proj_components := (x) -> piecewise((x >= ell_min) and (x <= ell_max), Project(components, RF, obj[parse("frame")]), 0);
  end if;

  if IsRod(obj) then
    if (proj_components(x)[2] <> 0) or (proj_components(x)[3] <> 0) then
      error "only axial loads are accepted in ROD objects"
    end if;
  end if;

  return table({
    parse("type")        = QFORCE,
    parse("components")  = proj_components,
    parse("coordinates") = [ell_min, ell_max],
    parse("target")      = obj[parse("name")]
    });
end proc: # MakeQForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsQForce := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a QFORCE object";

  if (obj[parse("type")] = QFORCE) and
     type(obj, table) and
     type(obj[parse("components")], procedure) and
     type(obj[parse("coordinates")], list) and
     type(obj[parse("target")], string) then
    return true;
  else
    return false;
  end if;
end proc: # IsQForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeQMoment := proc(
  components::{procedure,list(algebraic)}, # Distributed load components
  obj::BEAM,                               # Target object
  RF::FRAME := ground,                     # Reference frame in which the moment is defined
  {
    ell_min::algebraic := 0,                    # Initial application point (axial coordinate)
    ell_max::algebraic := obj[parse("length")]  # Final application point (axial coordinate)
  }, $)::QMOMENT;

  description "Define a QMOMENT object with inputs: distributed torque components, "
    "target object, initial and final application points (axial coordinates), "
    "optional reference frame in which the load components are defined "
    "(default = ground)";

  local proj_components, x;

  if type(components, procedure) then
    proj_components := unapply(Project(components(x), RF, obj[parse("frame")]),x);
  else
    proj_components := (x) -> piecewise((x >= ell_min) and (x <= ell_max), Project(components, RF, obj[parse("frame")]), 0);
  end if;

  return table({
    parse("type")        = QMOMENT,
    parse("components")  = proj_components,
    parse("coordinates") = [ell_min, ell_max],
    parse("target")      = obj[parse("name")]
    });
end proc: # MakeQMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsQMoment := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a QMOMENT object";

  if (obj[parse("type")] = QMOMENT) and
     type(obj, table) and
     type(obj[parse("components")], procedure) and
     type(obj[parse("coordinates")], list) and
     type(obj[parse("target")], string) then
    return true;
  else
    return false;
  end if;
end proc: # IsQMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeSupport := proc(
  name::string,            # Support name
  constrained_dof::list,   # Constrained degree of freedom
  objs::list({BEAM, ROD}), # Target objects
  ells::list(algebraic),   # Support locations
  RF::FRAME := ground,     # Reference frame of the support
  {
    stiffness::{procedure,list(algebraic)} := [ # Stiffness components (default = infinite)
      infinity, infinity, infinity,
      infinity, infinity, infinity
    ] *~ constrained_dof
  }, $)::SUPPORT;

  description "Define a SUPPORT object with inputs: support name, constrained "
    "degree of freedom, target objects, list of support locations, optional "
    "reference frame in which the support is defined (default = ground)";

  local S, J_tmp, i, j, sr_F_names, sr_F_values_tmp, sr_M_names, sr_M_values_tmp, S_stiffness, x;

  for i from 1 to nops(objs) do
    if IsRod(objs[i]) and (ells[i] <> 0) and (ells[i] <> objs[i][parse("length")]) then
      error "SUPPORT objects can only be applied at extremes of ROD objects"
    end if;
    if IsRod(objs[i]) and (constrained_dof[4..6] <> [0, 0, 0]) then
      error "ROD objects supports can only have translational constraints"
    end if;
  end do;

  if type(stiffness, procedure) then
    S_stiffness := unapply(stiffness(x) *~ constrained_dof, x);
    if (S_stiffness(x) <> stiffness(x)) then
      WARNING("stiffness components not corresponding to constrained_dof are ignored");
    end if;
  else
    S_stiffness := (x) -> stiffness *~ constrained_dof;
    if (S_stiffness(x) <> stiffness) then
      WARNING("stiffness components not corresponding to constrained_dof are ignored");
    end if;
  end if;

  S := table({
    parse("type")                     = SUPPORT,
    parse("constrained_dof")          = constrained_dof,
    parse("coordinates")              = [0, op(ells)],
    parse("name")                     = name,
    parse("frame")                    = RF,
    parse("targets")                  = [EARTH[parse("name")], op(GetNames(objs))],
    parse("variables")                = [],
    parse("forces")                   = [],
    parse("moments")                  = [],
    parse("constraint_loads")         = [],
    parse("support_reactions")        = [], # Expressed in support reference frame
    parse("stiffness")                = S_stiffness,
    parse("displacements")            = []
    });

  # Build the temporary joint
  J_tmp := MakeJoint(name, constrained_dof, [EARTH, op(objs)], S[parse("coordinates")], RF);

  S[parse("variables")]                := J_tmp[parse("variables")];
  S[parse("forces")]                   := J_tmp[parse("forces")];
  S[parse("moments")]                  := J_tmp[parse("moments")];
  S[parse("constraint_loads")]         := J_tmp[parse("constraint_loads")];

  # Retrieve support force reactions
  sr_F_names := [FX, FY, FZ];
  for i from 1 to nops(S[parse("forces")]) do
    if (S[parse("forces")][i][parse("target")] = EARTH[parse("name")]) then
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
    if (S[parse("moments")][i][parse("target")] = EARTH[parse("name")]) then
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

  return op(S);
end proc: # MakeSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsSupport := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a SUPPORT object";

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
    return true;
  else
    return false;
  end if;

end proc: # IsSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsCompliantSupport := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a SUPPORT object with compliant "
    "constraints";

  local found, i;

  if not IsSupport(obj) then
    error "object is not a SUPPORT";
    return false;
  end if;

  found := false;
  for i from 1 to nops(obj[parse("stiffness")]) do
    if (obj[parse("stiffness")][i] <> infinity) and (obj[parse("constrained_dof")][i] = 1) then
      found := true;
      break;
    end if;
  end do;

  return found;
end proc; # IsCompliantSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanSupport := proc(
  obj::SUPPORT, # Support to be cleaned
  $)::SUPPORT;

  description "Clean SUPPORT object <obj> internal variables";

  obj[parse("constraint_loads")]         := [];
  obj[parse("support_reactions")]        := [];
  obj[parse("displacements")]            := [];
  return obj;
end proc: # CleanSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeJoint := proc(
  name::string,                                   # Joint name
  constrained_dof::list,                          # Constrained degree of freedom
  objs::list({BEAM, ROD, SUPPORT, JOINT, EARTH}), # Target objects
  ells::list(algebraic),                          # Joint locations
  RF::FRAME := ground,                            # Reference frame
  $)::JOINT;

  description "Make a JOINT object with inputs: joint name, constrained "
    "degrees of freedom, target objects, joint locations, and optional "
    "reference frame in which the joint is defined (default = ground)";

  local J, i, jf_comp, jm_comp, jf_comp_obj, jm_comp_obj, jm_indets, jf_indets,
    constraint;

  for i from 1 to nops(objs) do
    if IsRod(objs[i]) and (ells[i] <> 0) and (ells[i] <> objs[i][parse("length")]) then
      error "JOINT objects can only be applied at extremes of ROD objects";
    end if;
    if IsRod(objs[i]) and (constrained_dof[4..6] <> [0, 0, 0]) then
      error "ROD objects supports can only have translational constraints";
    end if;
  end do;

    J := table({
      parse("type")                     = JOINT,
      parse("constrained_dof")          = constrained_dof,
      parse("coordinates")              = ells,
      parse("name")                     = name,
      parse("frame")                    = RF,
      parse("targets")                  = GetNames(objs),
      parse("variables")                = [],
      parse("forces")                   = [],
      parse("moments")                  = [],
      parse("constraint_loads")         = []
      });

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
    # Use the non admissible loads to build the loads constraint
    constraint := convert(
      Project(jf_comp, RF, objs[i][parse("frame")])
      .~ <op((-1*objs[i][parse("admissible_loads")][1..3]) +~ 1)>,
      list);
    constraint := remove(x -> x = 0, constraint);
    J[parse("constraint_loads")] := [
      op(J[parse("constraint_loads")]),
      op(constraint)
      ];
    # Extract the survived components
    jf_indets := indets(jf_comp);
    # Check if there are reactions
    if (jf_comp_obj <> [0, 0, 0]) then
      # Create the reaction force between joint and obj
      JF_||(name)||_||(objs[i][parse("name")]) := MakeForce(jf_comp_obj, ells[i], objs[i], objs[i][parse("frame")]);
      JF_||(objs[i][parse("name")])||_||(name) := MakeForce(-jf_comp_obj, 0, J, objs[i][parse("frame")]);
      # Update the output joint
      J[parse("variables")] := [
        op(J[parse("variables")]),
        op(jf_indets)
        ];
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
    if (jm_comp_obj <> [0, 0, 0]) then
      # Create the reaction force between joint and obj
      JM_||(name)||_||(objs[i][parse("name")]) := MakeMoment(jm_comp_obj, ells[i], objs[i], objs[i][parse("frame")]);
      JM_||(objs[i][parse("name")])||_||(name) := MakeMoment(-jm_comp_obj, 0, J, objs[i][parse("frame")]);
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

  return op(J);
end proc: # MakeJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsJoint := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the objecy <obj> is a JOINT object";

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
    return true;
  else
    return false;
  end if;
end proc: # IsJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanJoint := proc(
  obj::JOINT, # Object to be cleaned
  $)::JOINT;

  description "Clean JOINT object <obj> internal variables";

  obj[parse("constraint_loads")]         := [];
  return obj;
end proc: # CleanJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeRod := proc(
  name::string,        # Object name
  ell::algebraic,      # Length (m)
  RF::FRAME := ground, # Reference frame
  {
    area::{algebraic, procedure} := 0,   # Section area (m^2)
    material::MATERIAL           := NULL # Material
  }, $)::ROD;

  description "Create a ROD object with inputs: object name, reference "
  "frame, length, and optional section area and material";

  local area_proc;

  if type(area, procedure) then
    area_proc := area;
  else
    area_proc := (x) -> piecewise((x >= 0) and (x <= ell), area, 0);
  end if;

  return table({
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
end proc: # MakeRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsRod := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a ROD object";

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
    return true;
  else
    return false;
  end if;
end proc: # IsRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanRod := proc(
  obj::ROD, # Object to be cleaned
  $)::ROD;

  description "Clean ROD object <obj> internal variables";

  obj[parse("internal_actions")] := [];
  obj[parse("displacements")]    := [];
  return obj;
end proc: # CleanRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeBeam := proc(
  name::string,        # Object name
  ell::algebraic,      # Length (m)
  RF::FRAME := ground, # Reference frame
  {
    area::{algebraic, procedure}                     := 0,      # Section area (m^2)
    shear_stiff_factor::{list(algebraic), procedure} := [0, 0], # Shear stiffness factor
    material::MATERIAL                               := NULL,   # Material object
    I_xx::{algebraic, procedure}                     := 0,      # Section x-axis inertia (m^4)
    I_yy::{algebraic, procedure}                     := 0,      # Section y-axis inertia (m^4)
    I_zz::{algebraic, procedure}                     := 0       # Section z-axis inertia (m^4)
  }, $)::BEAM;

  description "Create a BEAM object with inputs: object name, reference "
    "frame, length, and optional section area, inertias on x-, y- and z-axis "
    "and material";

  local area_proc, shear_stiff_factor_proc, I_xx_proc, I_yy_proc, I_zz_proc;

  if type(area, procedure) then
    area_proc := area;
  else
    area_proc := (x) -> piecewise((x >= 0) and (x <= ell), area, 0);
  end if;

  if type(shear_stiff_factor, procedure) then
    shear_stiff_factor_proc := shear_stiff_factor;
  else
    shear_stiff_factor_proc := (x) -> piecewise((x >= 0) and (x <= ell), shear_stiff_factor, [0, 0]);
  end if;

  if type(I_xx, procedure) then
    I_xx_proc := I_xx;
  else
    I_xx_proc := (x) -> piecewise((x >= 0) and (x <= ell), I_xx, 0);
  end if;

  if type(I_yy, procedure) then
    I_yy_proc := I_yy;
  else
    I_yy_proc := (x) -> piecewise((x >= 0) and (x <= ell), I_yy, 0);
  end if;

  if type(I_zz, procedure) then
    I_zz_proc := I_zz;
  else
    I_zz_proc := (x) -> piecewise((x >= 0) and (x <= ell), I_zz, 0);
  end if;

  return table({
    parse("type")               = BEAM,
    parse("name")               = name,
    parse("length")             = ell,
    parse("area")               = area_proc,
    parse("shear_stiff_factor") = shear_stiff_factor_proc,
    parse("material")           = material,
    parse("inertias")           = [I_xx_proc, I_yy_proc, I_zz_proc],
    parse("frame")              = RF,
    parse("admissible_loads")   = [1, 1, 1, 1, 1, 1],
    parse("internal_actions")   = [],
    parse("displacements")      = []
    });
end proc: # MakeBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsBeam := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a BEAM object";

  if (obj[parse("type")] = BEAM) and
     type(obj, table) and
     type(obj[parse("name")], string) and
     type(obj[parse("length")], algebraic) and
     type(obj[parse("area")], procedure) and
     type(obj[parse("shear_stiff_factor")], procedure) and
     type(obj[parse("material")], MATERIAL) and
     type(obj[parse("inertias")], list(procedure)) and
     type(obj[parse("frame")], FRAME) and
     type(obj[parse("admissible_loads")], list) and
     type(obj[parse("internal_actions")], list) and
     type(obj[parse("displacements")], list) then
    return true;
  else
    return false;
  end if;
end proc: # IsBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanBeam := proc(
  obj::BEAM, # Object to be cleaned
  $)::BEAM;

  description "Clean BEAM object <obj> internal variables";

  obj[parse("internal_actions")] := [];
  obj[parse("displacements")]    := [];
  return obj;
end proc: # CleanBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeSpringDisplacement := proc(
  spring_load::algebraic, # load on the spring
  stiffness::procedure,   # spring stiffness
  $)::algebraic;

  description "Compute the displacement of a spring with stiffness "
    "<stiffness> and load <load>";

  local disp, x;

  disp := RealDomain[solve](spring_load = integrate(stiffness(x), x = 0..Dx), Dx);

  return disp;
end proc: # ComputeSpringDisplacement

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeSpringEnergy := proc(
  spring_load::algebraic, # load on the spring
  stiffness::procedure,   # spring stiffness
  $)::algebraic;

  description "Compute the potential energy of a spring with stiffness "
    "<stiffness> and subject to a load <load>";

  local P, disp, x;

  disp := ComputeSpringDisplacement(spring_load, stiffness);
  P := integrate(integrate(stiffness(x), x = 0..disp), x = 0..disp);
  
  return P;
end proc: # ComputeSpringEnergy

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeSupportInducedDisplacements := proc(
  obj::SUPPORT, # Support object
  exts::{ # External actions
    list({FORCE, MOMENT}),
    set( {FORCE, MOMENT})
  },
  $)::list;

  description "Compute the displacements of a support object <obj> induced "
    "from the external actions <exts>";

  local disp, loads_f, loads_m, loads, i, x;

  # Initialize load vectors
  loads_f := [0, 0, 0];
  loads_m := [0, 0, 0];
  # Initialize displacement vector
  disp := [0, 0, 0, 0, 0, 0];

  # Compute the load vector
  for i from 1 to nops(exts) do
    if (exts[i][parse("target")] <> obj[parse("name")]) then
      error "load target is not the support object";
    end if;
    if (exts[i][parse("type")] = FORCE) then
      loads_f := loads_f +~ exts[i][parse("components")];
    elif (exts[i][parse("type")] = MOMENT) then
      loads_m := loads_m +~ exts[i][parse("components")];
    end if;
  end do;
  loads := loads_f union loads_m;

  # Compute the displacements
  for i from 1 to 6 do
    disp[i] := ComputeSpringDisplacement(loads[i], (x -> obj[parse("stiffness")](x)[i]));
  end do;

  return disp;
end proc; # ComputeSupportInducedDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeSupportDisplacements := proc(
  obj::SUPPORT, # Support object
  $)::list;

  description "Compute the displacements of the support <obj> from its "
    "support reactions";

  local sup_disp, disp_vec, i, disp, x;

  # Increase printf indentation
  print_indent := print_indent + print_increment;

  sup_disp := [];
  disp_vec := [tx, ty, tz, rx, ry, rz];

  for i from 1 to 6 do 
    if (obj[parse("constrained_dofs")][i] = 1) then
      disp := ComputeSpringDisplacement(obj[parse("support_reactions")][i],
        (x -> obj[parse("stiffness")](x)[i]));
      sup_disp := sup_disp union [disp_vec[i] = disp];
    end if;
  end do;

  printf(
    "%*sMessage (in ComputeSupportDisplacements) updating %s %s's displacements...\n",
    print_indent, "", obj[parse("type")], obj[parse("name")]
    );
  obj[parse("displacements")] := sup_disp;
  printf("%*sDONE\n", print_indent, "");

  # Decrease printf indentation
  print_indent := print_indent - print_increment;

  return ``;
end proc: # ComputeSupportDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeStructure := proc(
  objs::{ # Structure objects
    list({BEAM, ROD, SUPPORT, JOINT}),
    set( {BEAM, ROD, SUPPORT, JOINT})
  },
  exts::{ # External actions
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  } := [],
  {
    hyper_vars::{list ,set} := [],
      # Hyperstatic variables
    hyper_disp::{list, set} := [seq(0, 1..nops(hyper_vars))],
      # Hyperstatic displacements
    dimensions::string := "3D",
      # Structure dimension ("2D" or "3D")
    verbose::boolean := false
      # Verbose mode
  }, $)::STRUCTURE;

  description "Create a STRUCTURE object with inputs: structure objects, "
    "external forces, moments or distributed loads, hyperstatic variables and "
    "displacements and structure dimension (""2D"" or ""3D"")";

  local num_dof, i, names, candidate_hyp_vars;

  # Increase printf indentation
  print_indent := print_indent + print_increment;

  # Check for duplicate names
  names := [];
  for i from 1 to nops(objs) do
    if member(objs[i][parse("name")], names) then
      error "duplicate names found on structure objects";
    end if;
    names := names union [objs[i][parse("name")]];
  end do;

  num_dof := ComputeDOF(objs, parse("dimensions") = dimensions, parse("verbose") = verbose);

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
      printf("%*sMessage (in MakeStructure) "
        "hyperstatic structure detected with %d overconstrained directions\n",
        print_indent, "", abs(num_dof));
    end if;
  elif (num_dof > 0) then
    error "not enough constraints in the structure";
  else
    printf("%*sMessage (in MakeStructure) isostatic structure detected\n", print_indent, "");
  end if;

  # Decrease printf indentation
  print_indent := print_indent - print_increment;

  return table({
    parse("type")                      = STRUCTURE,
    parse("objects")                   = objs,
    parse("external_actions")          = exts,
    parse("dof")                       = num_dof,
    parse("hyperstatic_variables")     = hyper_vars,
    parse("hyperstatic_displacements") = hyper_disp,
    parse("dimensions")                = dimensions,
    parse("support_reactions_solved")  = false,
    parse("internal_actions_solved")   = false,
    parse("displacement_solved")       = false
    });
end proc: # MakeStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsStructure := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a STRUCTURE object";

  if (obj[parse("type")] = STRUCTURE) then
    return true;
  else
    return false;
  end if;
end proc: # IsStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanStructure := proc(
  obj::STRUCTURE, # Object to be cleaned
  $)::STRUCTURE;

  description "Clean STRUCTURE object <obj> internal variables";

  local i;

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
end proc: # CleanStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeDOF := proc(
  objs::{ # Structure objects
    list({BEAM, ROD, SUPPORT, JOINT}),
    set( {BEAM, ROD, SUPPORT, JOINT})
  },
  {
    dimensions::string := "3D", # Structure dimension ("2D" or "3D")
    verbose::boolean   := false # Verbose mode
  }, $)::integer;

  description "Compute the degree of freedom of the input structure objects";

  local dof, objs_tmp, i, j, k, vertex, G;

  # Increase printf indentation
  print_indent := print_indent + print_increment;

  dof      := 0;
  objs_tmp := objs union [EARTH];

  # Built connections graph
  vertex := [];
  printf("%*sMessage (in ComputeDOF) checking structure connections...\n", print_indent, "");
  for i from 1 to nops(objs_tmp) do
    vertex := vertex union [objs_tmp[i][parse("name")]];
    end do;
  G := GraphTheory[Graph](vertex);
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
    printf("%*sDONE\n", print_indent, "");
  else
    error "unconnected elements detected in the structure";
  end if;

printf("%*sMessage (in ComputeDOF) computing degrees of freedom...\n", print_indent, "");
for i from 1 to nops(objs_tmp) do
  if IsBeam(objs_tmp[i]) then
      if (dimensions = "2D") then
        dof := dof + 3;
      else
        dof := dof + 6;
      end if;
    elif IsRod(objs_tmp[i]) then
      if (dimensions = "2D") then
        dof := dof + 3;
      else
        dof := dof + 5;
      end if;
    elif IsJoint(objs_tmp[i]) then
      if (dimensions = "2D") then
        dof := dof - add(objs_tmp[i][constrained_dof][k], k = [1,2,6]) * (nops(objs_tmp[i][parse("targets")]) - 1);
      else
        dof := dof - add(objs_tmp[i][constrained_dof][k], k = 1..6) * (nops(objs_tmp[i][parse("targets")]) - 1);
      end if;
    elif IsSupport(objs_tmp[i]) then
      if (dimensions = "2D") then
        dof := dof - add(objs_tmp[i][constrained_dof][k], k = [1,2,6]) * (nops(objs_tmp[i][parse("targets")]) - 1);
      else
        dof := dof - add(objs_tmp[i][constrained_dof][k], k = 1..6) * (nops(objs_tmp[i][parse("targets")]) - 1);
      end if;
    end if;
  end do;
  printf("%*sDONE\n", print_indent, "");
  printf("%*sMessage (in ComputeDOF) display degrees of freedom... DOF = %d\n",print_indent, "", dof);

  # Display graph
  if (verbose) then
    printf("%*sMessage (in ComputeDOF) display connections graph...\n", print_indent, "");
    print(GraphTheory[DrawGraph](G));
  end if;

  # Decrease printf indentation
  print_indent := print_indent - print_increment;
  
  return dof;
end proc: # ComputeDOF

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

NewtonEuler := proc(
  exts::{ # External actions
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  },
  obj::{BEAM, ROD, SUPPORT, JOINT}, # Object to compute the equilibrium
  pole,                             # Pole to compute the equilibrium # TODO: which type is this?
  {
    dimensions::string := "3D",                   # Structure dimension ("2D" or "3D")
    lim::algebraic        := obj[parse("length")] # Upper limit of the integration
  }, $)

  description "Compute the Newton-Euler static equilibrium equations given a "
    "set of external actions and the axial coordinate of the pole";

  local eq_T, eq_R, i, x;

  # Increase printf indentation
  print_indent := print_indent + print_increment;

  # 2D case
  if (dimensions = "2D") then

    eq_T := [0, 0];
    for i from 1 to nops(exts) do
      if (exts[i][parse("target")] = obj[parse("name")]) then
        if IsForce(exts[i]) then
          eq_T := eq_T + exts[i][parse("components")][1..2];
        elif IsQForce(exts[i]) then
          eq_T := eq_T + map(integrate, exts[i][parse("components")](x)[1..2], x = 0..lim);
        end if;
      else
        WARNING("Message (in NewtonEuler) %1 is not applied to %2", exts[i], obj);
      end if;
    end do;

    eq_R := [0];
    for i from 1 to nops(exts) do
      if exts[i][parse("target")] = obj[parse("name")] then
        if IsMoment(exts[i]) then
          eq_R := eq_R + exts[i][parse("components")][3];
        elif IsForce(exts[i]) then
          eq_R := eq_R + [exts[i][parse("components")][2]] *~ (exts[i][parse("coordinate")] - pole);
        elif IsQForce(exts[i]) then
          eq_R := eq_R + map(integrate, [exts[i][parse("components")](x)[2]*~(x-pole)], x = 0..lim);
        elif IsQMoment(exts[i]) then
          eq_R := eq_R + map(integrate, exts[i][parse("components")](x)[3], x = 0..lim);
        end if;
      else
        WARNING("Message (in NewtonEuler) %1 is not applied to %2", exts[i], obj);
      end if;
    end do;

  # 3D case
  elif (dimensions = "3D") then

    eq_T := [0, 0, 0];
    for i from 1 to nops(exts) do
      if exts[i][parse("target")] = obj[parse("name")] then
        if IsForce(exts[i]) then
          eq_T := eq_T + exts[i][parse("components")];
        elif IsQForce(exts[i]) then
          eq_T := eq_T + map(integrate, exts[i][parse("components")](x), x = 0..lim);
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
          eq_R := eq_R + [0, -exts[i][parse("components")][3], exts[i][parse("components")][2]]
            *~ (exts[i][parse("coordinate")] - pole);
        elif IsQForce(exts[i]) then
          eq_R := eq_R + map(integrate,
            [0, -exts[i][parse("components")](x)[3]*~(x-pole), exts[i][parse("components")](x)[2]*~(x-pole)],
            x = 0..lim);
        elif IsQMoment(FMQ[i]) then
          eq_R := eq_R + map(integrate,
            exts[i][parse("components")](x), x=0..lim);
        end if;
      else
        WARNING("Message (in NewtonEuler) %1 is not applied to %2", exts[i], obj);
      end if;
    end do;

  else
    error("invalid dimension detected");
  end if;

  # Decrease printf indentation
  print_indent := print_indent - print_increment;

  return eq_T union eq_R;
end proc: # NewtonEuler

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SolveStructure := proc(
  struct::STRUCTURE, # Structure object
  {
    compute_intact::boolean     := false, # Internal actions computation flag
    compute_disp::boolean       := false, # Displacement computation flag
    shear_contribution::boolean := false, # Shear contribution flag
    verbose::boolean            := false  # Verbose mode
  }, $)

    description "Solve the static equilibrium of a structure with inputs: the "
      "structure object, the compute internal action enabling flag, the compute "
      "displacement enabling flag, the shear contribution enabling flag and the "
      "verbose mode";

    local i, g_load, S_obj, S_ext, S_support, S_joint, S_con_forces, vars, sol,
      obj, x;

  # Increase printf indentation
  print_indent := print_indent + print_increment;

    # Parsing inputs
    S_obj        := {};
    S_ext        := {};
    S_support    := {};
    S_joint      := {};
    S_con_forces := {};
    vars         := [];
    for i from 1 to nops(struct[parse("objects")]) do
      obj := struct[parse("objects")][i];
      if IsBeam(obj) or IsRod(obj) then
        S_obj := S_obj union {obj};
        end if;
      if IsSupport(obj) then
        S_support    := S_support union {obj};
        S_con_forces := S_con_forces union obj[parse("forces")] union obj[parse("moments")];
        vars         := vars union obj[parse("variables")];
        end if;
      if IsJoint(obj) then
        S_joint      := S_joint union {obj};
        S_con_forces := S_con_forces union obj[parse("forces")] union obj[parse("moments")];
        vars         := vars union obj[parse("variables")];
        end if;
        unassign('obj');
    end do;

    S_ext := struct[parse("external_actions")];

    # Add gravity distributed load
    if (_gravity <> [0, 0, 0]) then
        for i from 1 to nops(S_obj) do
          if IsBeam(S_obj[i]) then
            g_load||(S_obj[i][parse("name")]) := MakeQForce(
              (x -> _gravity *~ S_obj[i][parse("area")](x) *~ S_obj[i][parse("material")][parse("density")]),
              S_obj[i],ground
              );
            S_ext := S_ext union {g_load||(S_obj[i][parse("name")])};
          end if;
        end do;
    end if;

    if (struct[dof] = 0) then
      # Solve isostatic structure
      printf("%*sMessage (in SolveStructure) solving the isostatic structure...\n", print_indent, "");
      sol := IsostaticSolver(
        S_obj union S_joint union S_support,
        S_ext union S_con_forces,
        vars, struct[parse("dimensions")], parse("verbose") = verbose
        );
      printf("%*sDONE\n", print_indent, "");
      if (verbose) then
        printf("Message (in SolveStructure) solutions:\n");
        print(<sol>);
      end if;
      # Update support reactions properties
      printf("%*sMessage (in SolveStructure) updating support reactions fields...\n", print_indent, "");
        for i from 1 to nops(S_support) do
          S_support[i][parse("support_reactions")] := [
            seq(lhs(S_support[i][parse("support_reactions")][j]) = subs(sol, rhs(S_support[i][parse("support_reactions")][j])),
            j = 1..nops(S_support[i][parse("support_reactions")]))
          ];
        end do;
      printf("%*sDONE\n", print_indent, "");
    elif (struct[dof] < 0) then
      # Solve hyperstatic structure
      if (nops(struct[parse("hyperstatic_variables")]) <> -struct[dof]) then
        error "mismatch in the structure degrees of freedom, check the hyper"
          "static variables of the structure and update the structure object";
      end if;
      printf("%*sMessage (in SolveStructure) solving the hyperstatic structure\n", print_indent, "");
      sol := HyperstaticSolver(
        S_obj union S_joint union S_support,
        S_ext union S_con_forces,
        vars,
        struct[parse("hyperstatic_variables")],
        struct[parse("hyperstatic_displacements")],
        parse("dimensions")         = struct[parse("dimensions")],
        parse("verbose")            = verbose,
        parse("shear_contribution") = shear_contribution
        );
      printf("%*sDONE\n", print_indent, "");
      if (verbose) then
        printf("%*sMessage (in SolveStructure) hyperstatic solver solution:\n", print_indent, "");
        print(<sol>);
      end if;
      # Update objects internal actions
      for obj in S_obj do
        obj[parse("internal_actions")] := subs(sol, obj[parse("internal_actions")]);
      end do;
      # Set internal actions computed flag
      struct[parse("internal_actions_solved")] := true; 
      # Update support reactions properties
      printf("%*sMessage (in SolveStructure) updating support reactions fields...\n", print_indent, "");
      for i from 1 to nops(S_support) do
      S_support[i][parse("support_reactions")] := [
        seq(lhs(S_support[i][parse("support_reactions")][j]) = subs(sol,rhs(S_support[i][parse("support_reactions")][j])),
        j = 1..nops(S_support[i][parse("support_reactions")]))
        ];
      end do;
      printf("%*sDONE\n", print_indent, "");
    end if;

  # Set support reactions solved flag
  struct[parse("support_reactions_solved")] := true;

  # Compute displacements
  if (compute_disp) and not struct[parse("displacement_solved")] then
    ComputeDisplacements(
      S_obj, S_ext union S_con_forces, sol, parse("dimensions") = struct[parse("dimensions")],
      parse("verbose") = verbose
      );
    # Set displacements computed flag
    struct[parse("displacement_solved")] := true;
  end if;

  # Compute internal actions
  if (compute_intact) and not struct[parse("internal_actions_solved")] then
    # FIXME: in case of Hyperstatic Structure, the internal actions are already computed
    ComputeInternalActions(
      S_obj, S_ext union S_con_forces, sol, parse("dimensions") = struct[parse("dimensions")],
      parse("verbose") = verbose
      );
    # Set internal actions computed flag
    struct[parse("internal_actions_solved")] := true;
  end if;

  # Decrease printf indentation
  print_indent := print_indent - print_increment;

  return struct;
end proc: # SolveStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

HyperstaticSolver := proc(
  objs::{ # Structural objects
    list({BEAM, ROD, SUPPORT, JOINT}),
    set( {BEAM, ROD, SUPPORT, JOINT})
  },
  exts::{ # External actions
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  },
  vars::list,       # Variables
  hyper_vars::list, # Hyperstatic variables
  hyper_disp::list, # Hyperstatic displacements
  {
    dimensions::string          := "3D",  # Dimensions ("2D" or "3D")
    shear_contribution::boolean := false, # Shear contribution
    verbose::boolean            := false  # Verbose mode
  }, $)

  description "Solve hyperstatic structure with inputs objects, external "
    "actions, variables, hyperstatic variables, hyperstatic displacements, "
    "dimensions and verbosity";

  local hyper_eq, hyper_load, hyper_comps, hyper_compliant_disp, hyper_support, i, obj,
        iso_vars, iso_sol, hyper_sol, sol, P, S_obj;

  # Increase printf indentation
  print_indent := print_indent + print_increment;

  # Parse input objects and find objects with internal actions property
    S_obj := {};
    for i from 1 to nops(objs) do
      if IsBeam(objs[i]) or IsRod(objs[i]) then
        S_obj := S_obj union {objs[i]};
      end if;
    end do;

  printf("%*sMessage (in HyperstaticSolver) solving the hyperstatic variables...\n", print_indent, "");
  # Create a solution as function of the hyperstatic variables
  iso_vars := [seq(
    `if`(member(vars[i], hyper_vars), NULL, vars[i]),
    i = 1..nops(vars))
    ];
  iso_sol := IsostaticSolver(objs, exts, iso_vars,  parse("dimensions") = dimensions, parse("verbose") = verbose);

  # Compute internal actions
  ComputeInternalActions(S_obj, exts, iso_sol, parse("dimensions") = dimensions, parse("verbose") = verbose);

  # Compute structure internal energy
  P := ComputePotentialEnergy(objs, parse("shear_contribution") = shear_contribution);

  hyper_eq := [];

  for i from 1 to nops(hyper_vars) do
    # Get the supports involved in the hyperstatic equation
    for obj in objs do
      if IsSupport(obj) and (member(hyper_vars[i], obj[parse("variables")])) then
        hyper_support := copy(obj);
        if (IsCompliantSupport(hyper_support)) then
          # Get the support loads related to the hyperstatic variable
          for hyper_load in hyper_support[parse("forces")] union hyper_support[parse("moments")] do
            if (hyper_load[parse("target")] =  hyper_support[parse("name")]) and 
                (has(hyper_load[parse("components")], hyper_vars[i])) then
              hyper_comps := eval(hyper_load[parse("components")] *~ map(has, hyper_load[parse("components")], hyper_vars[i]), [true=1,false=0]);
              # Create temporary load keeping only the components 
              # related to the hyperstatic variable
              if (IsForce(hyper_load)) then
                hyper_load := MakeForce(hyper_comps,0,hyper_support,hyper_support[parse("frame")]);
              elif (IsMoment(hyper_load)) then
                hyper_load := MakeMoment(hyper_comps,0,hyper_support,hyper_support[parse("frame")]);
              end if;
              break;
            end if;
          end do; 
          # Compute support induced displacements
          hyper_compliant_disp := Norm2(ComputeSupportInducedDisplacements(hyper_support, {hyper_load}));
          # Compose the hyperstatic equation, case compliant support
          hyper_eq := hyper_eq union [diff(P, hyper_vars[i]) = hyper_disp[i] + hyper_compliant_disp];
        else
          # Compose the hyperstatic equation, case rigid support
          hyper_eq := hyper_eq union [diff(P, hyper_vars[i]) = hyper_disp[i]];
        end if;
      end if;
    end do;
  end do;

  # Solve hyperstatic equations
  hyper_sol := op(solve(hyper_eq, hyper_vars));
  printf("%*sDONE\n", print_indent, "");

  # Decrease printf indentation
  print_indent := print_indent - print_increment;

  return hyper_sol union subs(hyper_sol, iso_sol);
end proc: # HyperstaticSolver

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputePotentialEnergy := proc(
  objs::{ # Structure objects
    list({BEAM, ROD, SUPPORT, JOINT}),
    set( {BEAM, ROD, SUPPORT, JOINT})
  },
  {
    shear_contribution := false # Add shear contribution to the potential energy
  }, $)

  description "Compute the internal potential energy of the structure";

  local obj, P, x;

  # Increase printf indentation
  print_indent := print_indent + print_increment;

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
      if shear_contribution then
        # Shear action Ty contribution
        if (member(Ty, map(lhs, obj[parse("internal_actions")]))) and
            (subs(obj[parse("internal_actions")](x), Ty(x)) <> 0) then
          P := P + integrate(
            subs(obj[parse("internal_actions")](x),
              obj[parse("shear_stiff_factor")](x)[1]*Ty(x)^2/(2*obj[parse("material")][parse("shear_modulus")]*obj[parse("area")](x))
            ), x = 0..obj[parse("length")]);
        end if;
        # Shear action Tz contribution
        if (member(Tz, map(lhs, obj[parse("internal_actions")]))) and
            (subs(obj[parse("internal_actions")](x), Tz(x)) <> 0) then
          P := P + integrate(
            subs(obj[parse("internal_actions")](x),
              obj[parse("shear_stiff_factor")](x)[2]*Tz(x)^2/(2*obj[parse("material")][parse("shear_modulus")]*obj[parse("area")](x))
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
    elif IsSupport(obj[i]) and IsCompliantSupport(obj[i]) then
      # Support reaction Fx contribution
      if (subs(obj[parse("support_reactions")], Fx) <> 0) and
          (obj[parse("stiffness")](x)[1] <> infinity) then
        P := P + subs(obj[parse("support_reactions")], ComputeSpringEnergy(Fx, (x -> obj[parse("stiffness")](x)[1])));
      end if;
      # Support reaction Fy contribution
      if (subs(obj[parse("support_reactions")], Fy) <> 0) and
          (obj[parse("stiffness")](x)[2] <> infinity) then
        P := P + subs(obj[parse("support_reactions")], ComputeSpringEnergy(Fy, (x -> obj[parse("stiffness")](x)[2])));
      end if;
      # Support reaction Fz contribution
      if (subs(obj[parse("support_reactions")], Fz) <> 0) and
          (obj[parse("stiffness")](x)[3] <> infinity) then
        P := P + subs(obj[parse("support_reactions")], ComputeSpringEnergy(Fz, (x -> obj[parse("stiffness")](x)[3])));
      end if;
      # Support reaction Mx contribution
      if (subs(obj[parse("support_reactions")], Mx) <> 0) and
          (obj[parse("stiffness")](x)[4] <> infinity) then
        P := P + subs(obj[parse("support_reactions")], ComputeSpringEnergy(Mx, (x -> obj[parse("stiffness")](x)[4])));
      end if;
      # Support reaction My contribution
      if (subs(obj[parse("support_reactions")], My) <> 0) and
          (obj[parse("stiffness")](x)[5] <> infinity) then
        P := P + subs(obj[parse("support_reactions")], ComputeSpringEnergy(My, (x -> obj[parse("stiffness")](x)[5])));
      end if;
      # Support reaction Mz contribution
      if (subs(obj[parse("support_reactions")], Mz) <> 0) and
          (obj[parse("stiffness")](x)[6] <> infinity) then
        P := P + subs(obj[parse("support_reactions")], ComputeSpringEnergy(Mz, (x -> obj[parse("stiffness")](x)[6])));
      end if;
    end if;
  end do;

  # Decrease printf indentation
  print_indent := print_indent - print_increment;

  return P;
end proc: # PotentialEnergy

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsostaticSolver := proc(
  objs::{ # Structural objects
    list({BEAM, ROD, SUPPORT, JOINT}),
    set( {BEAM, ROD, SUPPORT, JOINT})
  },
  exts::{ # External actions
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  },
  vars::list, # Variables to solve
  {
    dimensions::string := "3D", # Dimension ("2D" or "3D")
    verbose::boolean   := false # Verbose mode
  }, $)

  description "Solve the isostatic structure equilibrium equation system given "
    "the structure objects, the external actions and the variables to solve";

  local eq, i, j, active_ext, sol, A, B, rank_eq, vars_tmp;

  # Increase printf indentation
  print_indent := print_indent + print_increment;

  # Compute structure equations
  printf("%*sMessage (in IsostaticSolver) computing the equilibrium equation for "
    "the isostatic structure...\n", print_indent, "");
  eq := [];
    for i from 1 to nops(objs) do
    active_ext := {};
    for j from 1 to nops(exts) do
      if (exts[j][parse("target")] = objs[i][parse("name")]) then
        active_ext := active_ext union {exts[j]};
            end if;
        end do;
    eq := eq union NewtonEuler(active_ext, objs[i], 0, parse("dimensions") = dimensions);
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
  printf("%*sDONE\n", print_indent, "");

  # Structure equations check
  A, B    := LinearAlgebra[GenerateMatrix](eq, vars_tmp);
  rank_eq := LinearAlgebra[Rank](A);

  if (verbose) then
    printf("%*sMessage (in IsostaticSolver) structure equilibrium equations:\n", print_indent, "");
        print(<op(eq)>);
    printf("%*sMessage (in IsostaticSolver) structure unknown variables:\n", print_indent, "");
    print(vars_tmp);
    end if;

  if (rank_eq <> nops(vars_tmp)) then
    error "inconsistent system of equation, got %1 independent equations and "
      "%2 variables, check structure supports and joints",
      rank_eq, nops(vars_tmp);
  end if;

  # Solve structure equations
  printf("%*sMessage (in IsostaticSolver) computing the structure reaction forces...\n", print_indent, "");
  sol := simplify(op(solve(eq, vars_tmp)));
  printf("%*sDONE\n", print_indent, "");

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
  {
    dimensions::string := "3D", # Dimension ("2D" or "3D")
    verbose::boolean   := false # Verbose mode
  }, $)

  description "Programmatic computation of internal actions for structure"
    "objects with given external actions and structure solution";

  local i, j, active_ext, subs_ext;

  # Increase printf indentation
  print_indent := print_indent + print_increment;

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
    InternalActions(objs[i], active_ext, parse("dimensions") = dimensions);
  end do;

  # Decrease printf indentation
  print_indent := print_indent - print_increment;

end proc: # ComputeInternalActions

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

InternalActions := proc(
  obj::{BEAM, ROD}, # Structure object
  exts::{ # External actions
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  },
  {
    dimensions::string := "3D" # Dimension ("2D" or "3D")
  }, $)

  description "Programmatic computation of internal actions for structure "
    "object with given external actions, it returns the internal actions as "
    "function of the axial variable 'x'";

  local i, ia, N_sol, Ty_sol, Tz_sol, Mx_sol, My_sol, Mz_sol, x;

  # Increase printf indentation
  print_indent := print_indent + print_increment;

  # 2D case
  if (dimensions = "2D") then

  N_sol  := 0;
  Ty_sol := 0;
  Mz_sol := 0;

  # Compute internal actions for concentrated loads as effect overlay
  for i from 1 to nops(exts) do
    if IsForce(exts[i]) then
      N_sol  := `simplify/piecewise`(N_sol  - piecewise(x >= exts[i][parse("coordinate")] and x <= obj[parse("length")], exts[i][parse("components")][1]), x);
      Ty_sol := `simplify/piecewise`(Ty_sol + piecewise(x >= exts[i][parse("coordinate")] and x <= obj[parse("length")], exts[i][parse("components")][2]), x);
      Mz_sol := `simplify/piecewise`(Mz_sol + integrate(piecewise(x >= exts[i][parse("coordinate")] and x <= obj[parse("length")], exts[i][parse("components")][2]), x = 0..x), x);
    elif IsMoment(exts[i]) then
      Mz_sol := `simplify/piecewise`(Mz_sol - piecewise(x >= exts[i][parse("coordinate")] and x <= obj[parse("length")], exts[i][parse("components")][3]), x);
    elif IsQForce(exts[i]) then
      N_sol  := `simplify/piecewise`(N_sol  - integrate(exts[i][parse("components")](x)[1], x = 0..x), x);
      Ty_sol := `simplify/piecewise`(Ty_sol + integrate(exts[i][parse("components")](x)[2], x = 0..x), x);
      Mz_sol := `simplify/piecewise`(Mz_sol + integrate(integrate(exts[i][parse("components")](x)[2], x = 0..x), x = 0..x), x);
    elif IsQMoment(FMQ[i]) then
      Mz_sol := `simplify/piecewise`(Mz_sol - integrate(exts[i][parse("components")](x)[3], x = 0..x), x);
    end if;
    end do;

  ia := [
    N = unapply(N_sol, x), Ty = unapply(Ty_sol, x), Mz = unapply(Mz_sol, x)
            ];

  # 3D case
  elif (dimensions = "3D") then

    N_sol  := 0;
    Ty_sol := 0;
    Tz_sol := 0;
    Mx_sol := 0;
    My_sol := 0;
    Mz_sol := 0;

    # Compute internal actions for concentrated loads as effect overlay
    for i from 1 to nops(exts) do
      if IsForce(exts[i]) then
        N_sol  := `simplify/piecewise`(N_sol  - piecewise(x >= exts[i][parse("coordinate")] and x <= obj[parse("length")], exts[i][parse("components")][1]), x);
        Ty_sol := `simplify/piecewise`(Ty_sol + piecewise(x >= exts[i][parse("coordinate")] and x <= obj[parse("length")], exts[i][parse("components")][2]), x);
        Tz_sol := `simplify/piecewise`(Tz_sol + piecewise(x >= exts[i][parse("coordinate")] and x <= obj[parse("length")], exts[i][parse("components")][3]), x);
        My_sol := `simplify/piecewise`(My_sol + integrate(piecewise(x >= exts[i][parse("coordinate")] and x <= obj[parse("length")], exts[i][parse("components")][3]), x = 0..x), x);
        Mz_sol := `simplify/piecewise`(Mz_sol + integrate(piecewise(x >= exts[i][parse("coordinate")] and x <= obj[parse("length")], exts[i][parse("components")][2]), x = 0..x), x);
      elif IsMoment(exts[i]) then
        Mx_sol := `simplify/piecewise`(Mx_sol - piecewise(x >= exts[i][parse("coordinate")] and x <= obj[parse("length")], exts[i][parse("components")][1]), x);
        My_sol := `simplify/piecewise`(My_sol + piecewise(x >= exts[i][parse("coordinate")] and x <= obj[parse("length")], exts[i][parse("components")][2]), x);
        Mz_sol := `simplify/piecewise`(Mz_sol - piecewise(x >= exts[i][parse("coordinate")] and x <= obj[parse("length")], exts[i][parse("components")][3]), x);
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

  else
    error("invalid dimension detected");
  end if;

  if IsRod(obj) then
    ia := [ia[1]];
  end if;

  printf(
    "%*sMessage (in InternalActions) updating %s %s's internal actions...\n",
    print_indent, "", obj[parse("type")], obj[parse("name")]
    );
  obj[parse("internal_actions")] := ia;
  printf("%*sDONE\n", print_indent, "");

  # Decrease printf indentation
  print_indent := print_indent - print_increment;

  return ``;
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
    dimensions::string          := "3D", # Dimension ("2D" or "3D")
    verbose::boolean            := false, # Verbose mode
    shear_contribution::boolean := false # Add shear contribution to the potential energy
  }, $)

  description "Compute the Structure displacements";

  local dummy_Fx, dummy_Fy, dummy_Fz, dummy_Mx, dummy_My, dummy_Mz, obj,
        obj_copy, dummy_loads, subs_null_dummy, P, x, dFx, dFy, dFz, dMx, dMy, 
        dMz, disp;

  # Increase printf indentation
  print_indent := print_indent + print_increment;

  # Cicle on the structure objects
  for obj in objs do
    # Create a copy of the object
    obj_copy := copy(obj);

    # Beam
    if IsBeam(obj_copy) then
      # Create dummy loads
      dummy_Fx := MakeForce([dFx,0,0], x, obj_copy, obj_copy[parse("frame")]);
      dummy_Fy := MakeForce([0,dFy,0], x, obj_copy, obj_copy[parse("frame")]);
      dummy_Fz := MakeForce([0,0,dFz], x, obj_copy, obj_copy[parse("frame")]);
      dummy_Mx := MakeMoment([dMx,0,0], x, obj_copy, obj_copy[parse("frame")]);
      dummy_My := MakeMoment([0,dMy,0], x, obj_copy, obj_copy[parse("frame")]);
      dummy_Mz := MakeMoment([0,0,dMz], x, obj_copy, obj_copy[parse("frame")]);

      dummy_loads := {dummy_Fx, dummy_Fy, dummy_Fz, dummy_Mx, dummy_My, dummy_Mz};

      # Compute internal actions of the object copy
      ComputeInternalActions({obj_copy}, exts union dummy_loads, sol, parse("dimensions") = dimensions, parse("verbose") = verbose);

      # Compute object potential energy
      P := ComputePotentialEnergy({obj_copy}, parse("shear_contribution") = shear_contribution);

      # null dummy loads substitution list
      subs_null_dummy := [dFx, dFy, dFz, dMx, dMy, dMz] =~ [0,0,0,0,0,0];

      # Compute displacements
      disp := [0,0,0,0,0,0];
      disp[1] := tx = unapply(subs(subs_null_dummy, diff(P, dFx)),x);
      disp[2] := ty = unapply(subs(subs_null_dummy, diff(P, dFy)),x);
      disp[3] := tz = unapply(subs(subs_null_dummy, diff(P, dFz)),x);
      disp[4] := rx = unapply(subs(subs_null_dummy, diff(P, dMx)),x);
      disp[5] := ry = unapply(subs(subs_null_dummy, diff(P, dMy)),x);
      disp[6] := rz = unapply(subs(subs_null_dummy, diff(P, dMz)),x);

      # Update object displacements
      obj[parse("displacements")] := disp;

    # Rod
    elif IsRod(obj_copy) then
      # Create dummy loads
      dummy_Fx := MakeForce([dFx,0,0], x, obj_copy, obj_copy[parse("frame")]);

      dummy_loads := {dummy_Fx};

      # Compute internal actions of the object copy
      ComputeInternalActions({obj_copy}, exts union dummy_loads, sol, parse("dimensions") = dimensions, parse("verbose") = verbose);

      # Compute object potential energy
      P := ComputePotentialEnergy({obj_copy}, parse("shear_contribution") = shear_contribution);

      # null dummy loads substitution list
      subs_null_dummy := [dFx] =~ [0];

      # Compute displacements
      obj[parse("displacements")] := [tx = unapply(subs(subs_null_dummy, diff(P, dFx)),x)];

    # Support
    elif IsSupport(obj_copy) then
      # Compute displacements
      ComputeSupportDisplacements(obj);
    end if;
  end do;

  # Decrease printf indentation
  print_indent := print_indent - print_increment;

  return ``;
end proc: # ComputeDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotBeam := proc(
  obj::BEAM, # Beam to be plot
  $)

  description "Plot a the SUPPORT object <obj>";

  local P1, P2;

  P1 := Origin(obj[parse("frame")]);
  P2 := Origin(obj[parse("frame")].Translate(0, 0, obj[parse("length")]));

  return plots:-display(
    plottools:-line(convert(P1[1..3], list), convert(P2[1..3], list)),
    linestyle = solid, color = "SteelBlue");
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotRod := proc(
  obj::ROD, # Rod to be plot
  $)

  description "Plot a the ROD object <obj>";

  local P1, P2;

  P1 := Origin(obj[parse("frame")]);
  P2 := Origin(obj[parse("frame")].Translate(0, 0, obj[parse("length")]));

  return plots:-display(
    plottools:-line(convert(P1[1..3], list), convert(P2[1..3], list)),
    linestyle = dash, color = "Niagara DarkOrchid");
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotJoint := proc(
  obj::JOINT,            # Joint to be plot
  dim::algebraic := 0.1, # Dimensions
  $)

  description "Plot a the JOINT object <obj>";

  local O;

  O := Origin(
    parse(obj[parse("targets")][1])[parse("frame")].
    Translate(0,0,obj[parse("coordinates")][1])
    );

  return plots:-display(
    plottools:-sphere(convert(O[1..3], list), dim),
    linestyle = solid, color = "SteelBlue");
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotSupport := proc(
  obj::SUPPORT,          # Joint to be plot
  dim::algebraic := 0.1, # Dimensions
  $)

  local O;

  O := Origin(
    parse(obj[parse("targets")][2])[parse("frame")].
    Translate(0, 0, obj[parse("coordinates")][2])
    );

  return plots:-display(
    plottools:-sphere(convert(O[1..3], list), dim),
    linestyle = solid, color = "Niagara DarkOrchid");
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PlotStructure := proc(
  obj::STRUCTURE, # Structure to be plot
  $)

  description "Plot a the STRUCTURE object <obj>";

  local plt, i;

  plt := []:
  for i from 1 to nops(obj[parse("objects")]) do
    if (obj[parse("objects")][i][parse("type")] = BEAM) then
      plt := [op(plt), PlotBeam(obj[parse("objects")][i])];
    elif (obj[parse("objects")][i][parse("type")] = ROD) then
      plt := [op(plt), PlotRod(obj[parse("objects")][i])];
    elif (obj[parse("objects")][i][parse("type")] = SUPPORT) then
      plt := [op(plt), PlotSupport(obj[parse("objects")][i])];
    elif (obj[parse("objects")][i][parse("type")] = JOINT) then
      plt := [op(plt), PlotJoint(obj[parse("objects")][i])];
    end if;
  end do;
  return plt;
end proc:

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module:

# That's all folks!
