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

export  IsEarth,
        Show,
        Rotate,
        Translate,
        Project,
        Union,
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
        IsStructural,
        IsLoad,
        IsConstraint;

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
        lib_base_path;

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
  TypeTools[AddType](LOAD, IsLoad);
  TypeTools[AddType](CONSTRAINT, IsConstraint);
  TypeTools[AddType](STRUCTURAL, IsStructural);
  #TypeTools[AddType](TRUSSMEOBJECT, IsTrussMeObject);

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

  EARTH := table({
    type             = EARTH,
    name             = "earth",
    length           = 0,
    frame            = ground,
    admissible_loads = [1, 1, 1, 1, 1, 1]
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

  # Protect Module exported Functions
  protect(
    'Show',
    'Rotate',
    'Translate',
    'Project',
    'Union',
    'InverseFrame',
    'IsFrame',
    'IsEarth',
    'Origin',
    'Uvec',
    'UvecX',
    'UvecY',
    'UvecZ',
    'MakeMaterial',
    'IsMaterial',
    'MakeBeam',
    'IsBeam',
    'MakeRod',
    'IsRod',
    'MakeJoint',
    'IsJoint',
    'MakeSupport',
    'IsSupport',
    'IsCompliantSupport',
    'MakeForce',
    'IsForce',
    'MakeMoment',
    'IsMoment',
    'MakeQForce',
    'IsQForce',
    'MakeQMoment',
    'IsQMoment',
    'MakeStructure',
    'IsStructure',
    'SolveStructure',
    'IsStructural',
    'IsLoad',
    'IsConstraint'
  );

end proc: # Protect

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
     (obj[type] = EARTH) and
     (obj[length] = 0) and
     (obj[frame] = ground) and
     (obj[admissible_loads] = [1, 1, 1, 1, 1, 1]) then
    return true;
  else
    return false;
  end if:
end proc: # IsEarth

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Union := proc(
  A::{list, set}, # Object A to be united
  B::{list, set},  # Object B to be united
  $)::{list, set};

  description "Extension of union operator to list objects <A> and <B>";

  if type(A, 'set') and type(B, 'set') then
    return {op(A), op(B)};
  else
    return [op(A), op(B)];
  end if:
end proc: # Union

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Show := proc(
  tab::table, # Table to be shown
  $)::nothing;

  description "Show the content of a table <tab>";

  print(tab = tab[type](op(op(tab))));
end proc: # Show

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

GetNames := proc(
  objs::{ # Structural elements
    list({MATERIAL, STRUCTURAL}),
    set( {MATERIAL, STRUCTURAL})
  }, $)::{list(string), set(string)};

  description "Get names of a list/set of objects <objs>";

  if type(objs, 'list') then
    return {seq(objs[i][name], i = 1..nops(objs))};
  else
    return [seq(objs[i][name], i = 1..nops(objs))];
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
  axis::symbol,  # Rotation axis
  angle::scalar, # Rotation angle (rad)
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
  x::scalar, # X-axis translation component
  y::scalar, # Y-axis translation component
  z::scalar, # Z-axis translation component
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
  axis::symbol, # Axis of the unit vector
  RF::FRAME,    # Reference frame
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
  RF::FRAME, # Reference frame
  $)::vector;

  description "Extract the x-axis unit vector of the reference frame <RF>";

  return <RF[1,1], RF[2,1], RF[3,1], 0>;
end proc: # UvecX

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

UvecY := proc(
  RF::FRAME, # Reference frame
  $)::vector;

  description "Extract the y-axis unit vector of the reference frame <RF>";

  return <RF[1,2], RF[2,2], RF[3,2], 0>;
end proc: # UvecY

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

UvecZ := proc(
  RF::FRAME, # Reference frame
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
    name::string            := "steel",   # Name of the material
    elastic_modulus::scalar := 210.0E+09, # Elastic modulus (Pa)
    poisson_modulus::scalar := 0.3,       # Poisson modulus (-)
    shear_modulus::scalar   := 77.0E+09,  # Shear modulus (Pa)
    density::scalar         := 7.4E+03    # Density (kg/m^3)
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
  obj, # Object to be checked
  $)::boolean;

  description "Check if the input object <obj> is a MATERIAL object";

  if (obj[type] = MATERIAL) and
     type(obj, table) and
     type(obj[name], string) and
     type(obj[elastic_modulus], scalar) and
     type(obj[poisson_modulus], scalar) and
     type(obj[shear_modulus], scalar) and
     type(obj[density], scalar) then
    return true;
  else
    return false;
  end if;
end proc: # IsMaterial

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeForce := proc(
  components::list,         # Force components
  ell::scalar,              # Application point (axial coordinate)
  obj::{STRUCTURAL, EARTH}, # Target object
  RF::FRAME := ground,      # Reference frame in which the force is defined
  $)::FORCE;

description "Define a FORCE object with inputs: force components, force "
    "application axial coordinate <ell> in [0,L], target object <obj>, optional "
    "reference frame <RF> in which the force is defined (default = ground)";

  local proj_components;

  if IsBeam(obj) or IsRod(obj) then
    if (evalf(ell) < 0) or (evalf(ell) > evalf(obj[length])) then
      error "force application point must be in [0,L] range";
    end if;
  end if;

  proj_components := Project([op(components), 0], RF, obj[frame])[1..3];
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
    parse("target")     = obj[name]
    });
end proc: # MakeForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsForce := proc(
  obj, # Object to be checked
  $)::boolean;

  description "Check if obj is a FORCE object";

  if (obj[type] = FORCE) and
     type(obj, table) and
     type(obj[components], list) and
     type(obj[coordinate], scalar) and
     type(obj[target], string) then
    return true;
  else
    return false;
  end if;
end proc: # IsForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeMoment := proc(
  components::list,                   # Moment components
  ell::scalar,                        # Application point (axial coordinate)
  obj::{BEAM, SUPPORT, JOINT, EARTH}, # Target object
  RF::FRAME := ground,                # Reference frame in which the moment is defined
  $)::MOMENT;

  description "Define a MOMENT object with inputs: moment components, "
    "moment application axial coordinate [0,L], target object, optional "
    "reference frame in which the moment is defined (default = ground)";

  local proj_components;

  # FIXME: consider the case of symbolic length or ell
  if IsBeam(obj) then
    if (evalf(ell) < 0) or (evalf(ell) > evalf(obj[length])) then
      error "moment application point must be in [0,L] range";
    end if;
  end if;

  proj_components := Project([op(components), 0], RF, obj[frame])[1..3];
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
    parse("target")     = obj[name]
    });
end proc: # MakeMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsMoment := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a MOMENT object";

  if (obj[type] = MOMENT) and
     type(obj, table) and
     type(obj[components], list) and
     type(obj[coordinate], scalar) and
     type(obj[target], string) then
    return true;
  else
    return false;
  end if;
end proc: # IsMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeQForce := proc(
  components::list,    # Distributed load components
  obj::{BEAM, ROD},    # Target object
  RF::FRAME := ground, # Reference frame in which the object is defined
  {
    ell_min::scalar := 0,          # Initial application point (axial coordinate)
    ell_max::scalar := obj[length] # Final application point (axial coordinate)
  }, $)::QFORCE;

  description "Define a 'QFORCE' object with inputs: distributed load components, "
    "target object, initial and final application points (axial coordinates), "
    "optional reference frame in which the load components are defined "
    "(default = ground)";

  local proj_components, q_x, q_y, q_z;

  proj_components := Project([op(components), 0], RF, obj[frame])[1..3];

  if IsRod(obj) then
    if (proj_components[2] <> 0) or (proj_components[3] <> 0) then
      error "only axial loads are accepted in ROD objects"
    end if;
  end if;

  q_x := (x) -> piecewise((x >= ell_min) and (x <= ell_max), proj_components[1], 0);
  q_y := (x) -> piecewise((x >= ell_min) and (x <= ell_max), proj_components[2], 0);
  q_z := (x) -> piecewise((x >= ell_min) and (x <= ell_max), proj_components[3], 0);

  return table({
    parse("type")        = QFORCE,
    parse("components")  = [q_x, q_y, q_z],
    parse("coordinates") = [ell_min, ell_max],
    parse("target")      = obj[name]
    });
end proc: # MakeQForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsQForce := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if obj is a QFORCE object";

  if (obj[type] = QFORCE) and
     type(obj, table) and
     type(obj[components], list) and
     type(obj[coordinates], list) and
     type(obj[target], string) then
    return true;
  else
    return false;
  end if;
end proc: # IsQForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeQMoment := proc(
  components::list,    # Distributed load components
  obj::BEAM,           # Target object
  RF::FRAME := ground, # Reference frame in which the moment is defined
  {
    ell_min::scalar := 0,           # Initial application point (axial coordinate)
    ell_max::scalar := obj[length]  # Final application point (axial coordinate)
  }, $)::QMOMENT;

  description "Define a QMOMENT object with inputs: distributed torque components, "
    "target object, initial and final application points (axial coordinates), "
    "optional reference frame in which the load components are defined "
    "(default = ground)";

  local proj_components, q_x, q_y, q_z;

  proj_components := Project([op(components), 0], RF, obj[frame])[1..3];

  q_x := (x) -> piecewise((x >= ell_min) and (x <= ell_max), proj_components[1], 0);
  q_y := (x) -> piecewise((x >= ell_min) and (x <= ell_max), proj_components[2], 0);
  q_z := (x) -> piecewise((x >= ell_min) and (x <= ell_max), proj_components[3], 0);

  return table({
    parse("type")        = QMOMENT,
    parse("components")  = [q_x, q_y, q_z],
    parse("coordinates") = [ell_min, ell_max],
    parse("target")      = obj[name]
    });
end proc: # MakeQMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsQMoment := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if obj is a QMOMENT object";

  if (obj[type] = QMOMENT) and
     type(obj, table) and
     type(obj[components], list) and
     type(obj[coordinates], list) and
     type(obj[target], string) then
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
  ells::list,              # Support locations
  RF::FRAME := ground,     # Reference frame of the support
  {
    stiffness::list := [ # Stiffness components [Ktx, Kty, Ktz, Krx, Kry, Krz]
      infinity, infinity, infinity,
      infinity, infinity, infinity
    ]
  }, $)::SUPPORT;

  description "Define a SUPPORT object with inputs: support name, constrained "
    "degree of freedom, target objects, list of support locations, optional "
    "reference frame in which the support is defined (default = ground)";

  local S, J_tmp, i, j, sr_F_names, sr_F_values_tmp, sr_M_names, sr_M_values_tmp;

  for i from 1 to nops(objs) do
    if IsRod(objs[i]) and (ells[i] <> 0) and (ells[i] <> objs[i][length]) then
      error "SUPPORT objects can only be applied at extremes of ROD objects"
    end if;
    if IsRod(objs[i]) and (constrained_dof[4..6] <> [0, 0, 0]) then
      error "ROD objects supports can only have translational constraints"
    end if;
  end do;

  S := table({
    parse("type")                     = SUPPORT,
    parse("constrained_dof")          = constrained_dof,
    parse("coordinates")              = [0, op(ells)],
    parse("name")                     = name,
    parse("frame")                    = RF,
    parse("targets")                  = [EARTH[name], op(GetNames(objs))],
    parse("variables")                = [],
    parse("forces")                   = [],
    parse("moments")                  = [],
    parse("constraint_loads")         = [],
    parse("constraint_displacements") = [],
    parse("support_reactions")        = [], # Expressed in support reference frame
    parse("stiffness")                = stiffness,
    parse("displacements")            = []
    });

  # Build the temporary joint
  J_tmp := MakeJoint(name, constrained_dof, [EARTH, op(objs)], S[coordinates], RF);

  S[variables]                := J_tmp[variables];
  S[forces]                   := J_tmp[forces];
  S[moments]                  := J_tmp[moments];
  S[constraint_loads]         := J_tmp[constraint_loads];
  S[constraint_displacements] := J_tmp[constraint_displacements];

  # Retrieve support force reactions
  sr_F_names := [FX, FY, FZ];
  for i from 1 to nops(S[forces]) do
    if (S[forces][i][target] = "earth") then
      # Project forces in the support reference frame
      sr_F_values_tmp := Project(S[forces][i][components], ground, S[frame]);
      for j from 1 to 3 do
        if (sr_F_values_tmp[j] <> 0) then
          S[support_reactions] := [
            op(S[support_reactions]),
            sr_F_names[j] = -sr_F_values_tmp[j]
            ];
        end if;
      end do;
      break;
    end if;
  end do;

  # Retrieve support moments reactions
  sr_M_names := [MX, MY, MZ];
  for i from 1 to nops(S[moments]) do
    if (S[moments][i][target] = "earth") then
      # Project moments in the support reference frame
      sr_M_values_tmp := Project(S[moments][i][components], ground, S[frame]);
      for j from 1 to 3 do
        if (sr_M_values_tmp[j] <> 0) then
          S[support_reactions] := [
            op(S[support_reactions]),
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

  if (obj[type] = SUPPORT) and
     type(obj[constrained_dof], list) and
     type(obj[coordinates], list) and
     type(obj[name], string) and
     type(obj[frame], FRAME) and
     type(obj[targets], list(string)) and
     type(obj[variables], list) and
     type(obj[forces], list) and
     type(obj[moments], list) and
     type(obj[constraint_loads], list) and
     type(obj[constraint_displacements], list) and
     type(obj[support_reactions], list) and
     type(obj[stiffness], list) then
    return true;
  else
    return false;
  end if;

end proc: # IsSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsCompliantSupport := proc(
  obj, # Object to be checked
  $)

  description "Check if obj is a SUPPORT object with compliant constraints";

  local found;

  if not IsSupport(obj) then
    error "Object is not a SUPPORT";
    return false;
  end if;

  found := false;
  for i from 1 to nops(obj[stiffness]) do
    if not (obj[stiffness][i] = infinity) then
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

  obj[constraint_displacements] := [];
  return obj;
end proc: # CleanSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeJoint := proc(
  name::string,            # Joint name
  constrained_dof::list,   # Constrained degree of freedom
  objs::list({BEAM, ROD}), # Target objects
  ells::list,              # Joint locations
  RF::FRAME := ground,     # Reference frame in which the joint is defined
  $)

  description "Make a JOINT object with inputs: joint name, constrained "
    "degrees of freedom, target objects, joint locations, and optional "
    "reference frame in which the joint is defined (default = ground)";

  local J, i, jf_comp, jm_comp, jf_comp_obj, jm_comp_obj, jm_indets, jf_indets,
    constraint;

  for i from 1 to nops(objs) do
    if IsRod(objs[i]) and (ells[i] <> 0) and (ells[i] <> objs[i][length]) then
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
      parse("constraint_loads")         = [],
      parse("constraint_displacements") = []
      });

  # Add all the bodies forces
  for i from 1 to nops(objs) do
    # Create force compatible with the joint constrained dof
    jf_comp := convert(<
      JFx_||(J[name])||_||(objs[i][name]),
      JFy_||(J[name])||_||(objs[i][name]),
      JFz_||(J[name])||_||(objs[i][name])
      > *~ <op(constrained_dof[1..3])>,
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
      JF_||(name)||_||(objs[i][name]) := MakeForce(jf_comp_obj, ells[i], objs[i], objs[i][frame]);
      JF_||(objs[i][name])||_||(name) := MakeForce(-jf_comp_obj, 0, J, objs[i][frame]);
      # Update the output joint
      J[variables] := [
        op(J[variables]),
        op(jf_indets)
        ];
      J[forces] := [
        op(J[forces]),
        JF_||(name)||_||(objs[i][name]),
        JF_||(objs[i][name])||_||(name)
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
      > *~ <op(constrained_dof[4..6])>,
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
      JM_||(name)||_||(objs[i][name]) := MakeMoment(jm_comp_obj, ells[i], objs[i], objs[i][frame]);
      JM_||(objs[i][name])||_||(name) := MakeMoment(-jm_comp_obj, 0, J, objs[i][frame]);
      # Update the output joint
      J[variables] := [
        op(J[variables]),
        op(jm_indets)
        ];
      J[moments] := [
        op(J[moments]),
        JM_||(name)||_||(objs[i][name]),
        JM_||(objs[i][name])||_||(name)
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

  if (obj[type] = JOINT) and
     type(obj, table) and
     type(obj[constrained_dof], list) and
     type(obj[coordinates], list) and
     type(obj[name], string) and
     type(obj[frame], FRAME) and
     type(obj[targets], list(string)) and
     type(obj[variables], list) and
     type(obj[forces], list) and
     type(obj[moments], list) and
     type(obj[constraint_loads], list) and
     type(obj[constraint_displacements], list) then
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

  obj[constraint_displacements] := [];
  return obj;
end proc: # CleanJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeRod := proc(
  name::string, # Object name
  RF::FRAME,    # Reference frame
  ell::scalar,  # Length (m)
  {
    area::scalar       := 0,   # Section area (m^2)
    material::MATERIAL := NULL # Material
  }, $)::ROD;

  description "Create a ROD object with inputs: object name, reference "
  "frame, length, and optional section area and material";

  return table({
    parse("type")             = ROD,
    parse("name")             = name,
    parse("length")           = ell,
    parse("area")             = area,
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

  description "Check if obj is a ROD object";

  if (obj[type] = ROD) and
     type(obj, table) and
     type(obj[name], string) and
     type(obj[length], scalar) and
     type(obj[area], scalar) and
     type(obj[material], MATERIAL) and
     type(obj[frame], FRAME) and
     type(obj[admissible_loads], list) and
     type(obj[internal_actions], list) then
    return true;
  else
    return false;
  end if;
end proc: # IsRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanRod := proc(
  obj::{ROD}, # Object to be cleaned
  $)

  description "Clean ROD object internal variables";

  obj[internal_actions] := [];
  return obj;
end proc: # CleanRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeBeam := proc(
  name::string, # Object name
  RF::FRAME,    # Reference frame
  ell::scalar,  # Length (m)
  {
    area::scalar       := 0,    # Section area (m^2)
    material::MATERIAL := NULL, # Material object
    I_xx::scalar       := 0,    # Section x-axis inertia (m^4)
    I_yy::scalar       := 0,    # Section y-axis inertia (m^4)
    I_zz::scalar       := 0     # Section z-axis inertia (m^4)
  }, $)::BEAM;

  description "Create a BEAM object with inputs: object name, reference "
    "frame, length, and optional section area, inertias on x-, y- and z-axis "
    "and material";

  return table({
    parse("type")             = BEAM,
    parse("name")             = name,
    parse("length")           = ell,
    parse("area")             = area,
    parse("material")         = material,
    parse("inertias")         = [I_xx, I_yy, I_zz],
    parse("frame")            = RF,
    parse("admissible_loads") = [1, 1, 1, 1, 1, 1],
    parse("internal_actions") = [],
    parse("displacements")    = []
    });
end proc: # MakeBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsBeam := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a BEAM object";

  if (obj[type] = BEAM) and
     type(obj, table) and
     type(obj[name], string) and
     type(obj[length], scalar) and
     type(obj[area], scalar) and
     type(obj[material], MATERIAL) and
     type(obj[inertias], list) and
     type(obj[frame], FRAME) and
     type(obj[admissible_loads], list) and
     type(obj[internal_actions], list) then
    return true;
  else
    return false;
  end if;
end proc: # IsBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanBeam := proc(
  obj::BEAM, # Object to be cleaned
  $)::BEAM;

  description "Clean BEAM object internal variables";

  obj[internal_actions] := [];
  return obj;
end proc: # CleanBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsStructural := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a BEAM, ROD, SUPPORT or JOINT "
    "object";

  if IsBeam(obj) or IsRod(obj) or IsSupport(obj) or IsJoint(obj) then
    return true;
  else
    return false;
  end if;
end proc: # IsStructural

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsConstraint := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a SUPPORT or JOINT object";

  if IsSupport(obj) or IsJoint(obj) then
    return true;
  else
    return false;
  end if;
end proc: # IsConstraint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsLoad := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if the object <obj> is a FORCE, MOMENT, QFORCE or QMOMENT "
    "object";

  if IsForce(obj) or IsMoment(obj) or IsQForce(obj) or IsQMoment(obj) then
    return true;
  else
    return false;
  end if;
end proc: # IsLoad

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MakeStructure := proc(
  objs::{list(STRUCTURAL), set(STRUCTURAL)}, # Structural objects
  ext::{list(LOAD), set(LOAD)} := [],        # Load objects
  {
    hyper_vars::{list ,set} := [],
      # Hyperstatic variables
    hyper_disp::{list, set} := [seq(0, 1..nops(hyper_vars))],
      # Hyperstatic displacements
    dim::{string} := "3D"
      # Structure dimension ("2D" or "3D")
  }, $)::STRUCTURE;

  description "Create a STRUCTURE object with inputs: structure objects, "
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

  num_dof := ComputeDOF(objs, parse("dim") = dim);

  if (num_dof < 0) then
    if (nops(hyper_vars) <> -num_dof) then
      candidate_hyp_vars := [];
      for i from 1 to nops(objs) do
        if IsSupport(objs[i]) or IsJoint(objs[i]) then
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
      "field in an already existing STRUCTURE object by choosing from the "
      "folloving hyperstatic candidate variables: %2",
      abs(num_dof), candidate_hyp_vars
      );
    else
      printf("Message (in MakeStructure) "
        "hyperstatic structure detected with %d overconstrained directions\n",
        abs(num_dof));
    end if;
  elif (num_dof > 0) then
    error "not enough constraints in the structure";
  else
    printf("Message (in MakeStructure) isostatic structure detected");
  end if;

  return table({
    parse("type")                      = STRUCTURE,
    parse("objects")                   = objs,
    parse("external_actions")          = ext,
    parse("dof")                       = num_dof,
    parse("hyperstatic_variables")     = hyper_vars,
    parse("hyperstatic_displacements") = hyper_disp,
    parse("dimensions")                = dim,
    parse("support_reactions_solved")  = false,
    parse("internal_actions_solved")   = false,
    parse("displacement_solved")       = false
    });
end proc: # MakeStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsStructure := proc(
  obj::anything, # Object to be checked
  $)::boolean;

  description "Check if obj is a STRUCTURE object";

  if (obj[type] = STRUCTURE) then
    return true;
  else
    return false;
  end if;
end proc: # IsStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CleanStructure := proc(
  obj::STRUCTURE, # Object to be cleaned
  $)::STRUCTURE;

  description "Clean STRUCTURE object internal variables";

  local i;

  for i from 1 to nops(obj[objects]) do
    if IsBeam(obj[i]) then
      obj[objects][i] := CleanBeam(i);
    elif IsRod(obj[i]) then
      obj[objects][i] := CleanRod(i);
    elif IsSupport(obj[i]) then
      obj[objects][i] := CleanSupport(i);
    elif IsJoint(obj[i]) then
      obj[objects][i] := CleanJoint(i);
    end if;
  end do;
end proc: # CleanStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeDOF := proc(
  objs::{list(STRUCTURAL), set(STRUCTURAL)}, # Structure objects
  {
    dim::{string} := "3D" # Structure dimension ("2D" or "3D")
  }, $)::integer;

  description "Compute the degree of freedom of the input structure objects";

  local dof, objs_tmp, i, j, k, vertex, G;

  dof      := 0;
  objs_tmp := [op(objs), EARTH];

  # Built connections graph
  vertex := [];
  printf("Message (in ComputeDOF) checking structure connections... ");
  for i from 1 to nops(objs_tmp) do
    vertex := [op(vertex), objs_tmp[i][name]];
    end do;
  G := GraphTheory[Graph](vertex);
  for i from 1 to nops(objs_tmp) do
    if IsSupport(objs_tmp[i]) or IsJoint(objs_tmp[i]) then
      for j from 1 to nops(objs_tmp) do
        if (member(objs_tmp[j][name], objs_tmp[i][targets])) then
          GraphTheory[AddEdge](G, {objs_tmp[i][name], objs_tmp[j][name]});
        end if;
      end do;
    end if;
  end do;

  # Check graph connections
  if GraphTheory[IsConnected](G) then
    printf("DONE\n");
  else
    error "unconnected elements detected in the structure";
  end if;

printf("Message (in ComputeDOF) computing degrees of freedom... ");
for i from 1 to nops(objs_tmp) do
  if IsBeam(objs_tmp[i]) then
      if (dim = "2D") then
        dof := dof + 3;
      else
        dof := dof + 6;
      end if;
    elif IsRod(objs_tmp[i]) then
      if (dim = "2D") then
        dof := dof + 3;
      else
        dof := dof + 5;
      end if;
    elif IsJoint(objs_tmp[i]) then
      if (dim = "2D") then
        dof := dof - add(objs_tmp[i][constrained_dof][k], k = [1,2,6]) * (nops(objs_tmp[i][targets]) - 1);
      else
        dof := dof - add(objs_tmp[i][constrained_dof][k], k = 1..6) * (nops(objs_tmp[i][targets]) - 1);
      end if;
    elif IsSupport(objs_tmp[i]) then
      if (dim = "2D") then
        dof := dof - add(objs_tmp[i][constrained_dof][k], k = [1,2,6]) * (nops(objs_tmp[i][targets]) - 1);
      else
        dof := dof - add(objs_tmp[i][constrained_dof][k], k = 1..6) * (nops(objs_tmp[i][targets]) - 1);
      end if;
    end if;
  end do;
  printf("DONE\n");
  printf("Message (in ComputeDOF) display degrees of freedom... DOF = %d\n", dof);

  # Display graph
  print(GraphTheory[DrawGraph](G));

  return dof;
end proc: # ComputeDOF

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

NewtonEuler := proc(
  ext::{list(LOAD), set(LOAD)}, # External actions
  obj::STRUCTURAL,              # Object to compute the equilibrium
  pole,                         # Pole to compute the equilibrium # TODO: which type is this?
  {
    dim::{string} := "3D", # Structure dimension ("2D" or "3D")
    lim := obj[length]     # Upper limit of the integration
  }, $)

  description "Compute the Newton-Euler static equilibrium equations given a "
    "set of external actions and the axial coordinate of the pole";

  local eq_T, eq_R, i;

    # 2D case
  if (dim = "2D") then

    eq_T := [0, 0];
    for i from 1 to nops(ext) do
      if (ext[i][target] = obj[name]) then
        if IsForce(ext[i]) then
          eq_T := eq_T + ext[i][components][1..2];
        elif IsQForce(ext[i]) then
          eq_T := eq_T + map(integrate, ext[i][components][1..2](x), x = 0..lim);
        end if;
      else
        WARNING("Message (in NewtonEuler) %1 is not applied to %2", ext[i], obj);
      end if;
    end do;

    eq_R := [0];
    for i from 1 to nops(ext) do
      if ext[i][target] = obj[name] then
        if IsMoment(ext[i]) then
          eq_R := eq_R + ext[i][components][3];
        elif IsForce(ext[i]) then
          eq_R := eq_R + [ext[i][components][2]] *~ (ext[i][coordinate] - pole);
        elif IsQForce(ext[i]) then
          eq_R := eq_R + map(integrate, [ext[i][components][2](x)*~(x-pole)], x = 0..lim);
        elif IsQMoment(ext[i]) then
          eq_R := eq_R + map(integrate, ext[i][components][3](x), x = 0..lim);
        end if;
      else
        WARNING("Message (in NewtonEuler) %1 is not applied to %2", ext[i], obj);
      end if;
    end do;

  # 3D case
  elif (dim = "3D") then

    eq_T := [0, 0, 0];
    for i from 1 to nops(ext) do
      if ext[i][target] = obj[name] then
        if IsForce(ext[i]) then
          eq_T := eq_T + ext[i][components];
        elif IsQForce(ext[i]) then
          eq_T := eq_T + map(integrate, ext[i][components](x), x = 0..lim);
        end if;
      else
        WARNING("Message (in NewtonEuler) %1 is not applied to %2", ext[i], obj);
      end if;
    end do;

    eq_R := [0, 0, 0];
    for i from 1 to nops(ext) do
      if (ext[i][target] = obj[name]) then
        if IsMoment(ext[i]) then
          eq_R := eq_R + ext[i][components];
        elif IsForce(ext[i]) then
          eq_R := eq_R + [0, -ext[i][components][3], ext[i][components][2]]
            *~ (ext[i][coordinate] - pole);
        elif IsQForce(ext[i]) then
          eq_R := eq_R + map(integrate,
            [0, -ext[i][components][3](x)*~(x-pole), ext[i][components][2](x)*~(x-pole)],
            x = 0..lim);
        elif IsQMoment(FMQ[i]) then
          eq_R := eq_R + map(integrate,
            ext[i][components](x), x=0..lim);
        end if;
      else
        WARNING("Message (in NewtonEuler) %1 is not applied to %2", ext[i], obj);
      end if;
    end do;

  else
    error("invalid dimension detected");
  end if;

  return [op(eq_T), op(eq_R)];
end proc: # NewtonEuler

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

SolveStructure := proc(
  struct::STRUCTURE, # Structure object
  {
    compute_intact::{boolean}     := false, # Internal actions computation flag
    compute_disp::{boolean}       := false, # Displacement computation flag
    shear_contribution::{boolean} := false, # Shear contribution flag
    verbose::{boolean}            := false  # Verbose mode
  }, $)

    description "Solve the static equilibrium of a structure with inputs: the "
      "structure object, the compute internal action enabling flag, the compute "
      "displacement enabling flag, the shear contribution enabling flag and the "
      "verbose mode";

    local i, g_load, S_obj, S_ext, S_support, S_joint, S_con_forces, vars, sol,
      obj;

    # Parsing inputs
    S_obj        := {};
    S_ext        := {};
    S_support    := {};
    S_joint      := {};
    S_con_forces := {};
    vars         := [];
    for i from 1 to nops(struct[objects]) do
      obj := struct[objects][i];
      if IsBeam(obj) or IsRod(obj) then
        S_obj := {op(S_obj), obj};
        end if;
      if IsSupport(obj) then
        S_support    := {op(S_support), obj};
        S_con_forces := {op(S_con_forces), op(obj[forces]), op(obj[moments])};
        vars         := [op(vars), op(obj[variables])];
        end if;
      if IsJoint(obj) then
        S_joint      := {op(S_joint), obj};
        S_con_forces := {op(S_con_forces), op(obj[forces]), op(obj[moments])};
        vars         := [op(vars), op(obj[variables])];
        end if;
        unassign('obj');
    end do;

    S_ext := struct[external_actions];

    # Add gravity distributed load
    if (_gravity <> [0, 0, 0]) then
        for i from 1 to nops(S_obj) do
          if IsBeam(S_obj[i]) then
            g_load||(S_obj[i][name]) := MakeQForce(
              _gravity *~ S_obj[i][area] *~ S_obj[i][material][density],
              S_obj[i],ground
              );
            S_ext := {op(S_eg_load||(S_obj[i][name]))};
          end if;
        end do;
    end if;

    if (struct[dof] = 0) then
      # Solve isostatic structure
      printf("Message (in SolveStructure) solving the isostatic structure...\n");
      sol := IsostaticSolver(
        {op(S_obj), op(S_joint), op(S_support)},
        {op(S_ext), op(S_con_forces)},
        vars, struct[dimensions], parse("verbose") = verbose
        );
      printf("DONE\n");
      if (verbose) then
        printf("Message (in SolveStructure) solutions:\n");
        print(<sol>);
      end if;
      # Update support reactions properties
      printf("Message (in SolveStructure) updating support reactions fields... ");
        for i from 1 to nops(S_support) do
          S_support[i][support_reactions] := [
            seq(lhs(S_support[i][support_reactions][j]) = subs(sol, rhs(S_support[i][support_reactions][j])),
            j = 1..nops(S_support[i][support_reactions]))
          ];
        end do;
      printf("DONE\n");
    elif (struct[dof] < 0) then
      # Solve hyperstatic structure
      if (nops(struct[hyperstatic_variables]) <> -struct[dof]) then
        error "mismatch in the structure degrees of freedom, check the hyper"
          "static variables of the structure and update the structure object";
      end if;
      printf("Message (in SolveStructure) solving the hyperstatic structure\n");
      sol := HyperstaticSolver(
        {op(S_obj), op(S_joint), op(S_support)},
        {op(S_ext), op(S_con_forces)},
        vars,
        struct[hyperstatic_variables],
        struct[hyperstatic_displacements],
        parse("dim")      = struct[dimensions],
        parse("verbose")  = verbose,
        shear_contrib = shear_contribution
        );
      printf("DONE\n");
      if (verbose) then
        printf("Message (in SolveStructure) hyperstatic solver solution:\n");
        print(<sol>);
      end if;
      # Set internal actions computed flag
      struct[internal_actions_computed] := true;
      # Update support reactions properties
      printf("Message (in SolveStructure) updating support reactions fields... ");
      for i from 1 to nops(S_support) do
      S_support[i][support_reactions] := [
        seq(lhs(S_support[i][support_reactions][j]) = subs(sol,rhs(S_support[i][support_reactions][j])),
        j = 1..nops(S_support[i][support_reactions]))
        ];
      end do;
      printf("DONE\n");
    end if;

  # Set support reactions solved flag
  struct[support_reactions_solved] := true;

  # Compute displacements
  if (compute_displacements) and not struct[displacements_computed] then
    ComputeDisplacements(
      S_obj, {op(S_ext), op(S_con_forces)}, sol, struct[dimensions],
      parse("verbose") = verbose
      );
    # Set displacements computed flag
    struct[displacements_computed] := true;
  end if;

  # Compute internal actions
  if (compute_intact) and not struct[internal_actions_computed] then
    # FIXME: in case of Hyperstatic Structure, the internal actions are already computed
    ComputeInternalActions(
      S_obj, {op(S_ext), op(S_con_forces)}, sol, struct[dimensions],
      parse("verbose") = verbose
      );
    # Set internal actions computed flag
    struct[internal_actions_computed] := true;
  end if;

  return struct;
end proc: # SolveStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

HyperstaticSolver := proc(
  objs::{list(STRUCTURAL), set( STRUCTURAL)}, # Structural objects
  ext::{list(LOAD), set(LOAD)},               # External actions
  vars::{list},                               # Variables
  hyper_vars::{list},                         # Hyperstatic variables
  hyper_disp::{list},                         # Hyperstatic displacements
  {
    dim::{string}                 := "3D",  # Dimensions ("2D" or "3D")
    shear_contribution::{boolean} := false, # Shear contribution
    verbose::{boolean}            := false  # Verbose mode
  },
  $)

  description "Solve hyperstatic structure with inputs objects, external "
    "actions, variables, hyperstatic variables, hyperstatic displacements, "
    "dimensions and verbosity";

  local hyper_eq, i, iso_vars, iso_sol, hyper_sol, sol, P, S_obj;

  # Parse input objects and find objects with internal actions property
    S_obj := {};
    for i from 1 to nops(objs) do
      if IsBeam(objs[i]) or IsRod(objs[i]) then
        S_obj := {op(S_obj), objs[i]};
      end if;
    end do;

  printf("Message (in HyperstaticSolver) solving the hyperstatic variables... ");
  # Create a solution as function of the hyperstatic variables
  iso_vars := [seq(
    `if`(member(vars[i], hyper_vars), NULL, vars[i]),
    i = 1..nops(vars))
    ];
  iso_sol := IsostaticSolver(objs, ext, iso_vars, dim, parse("verbose") = verbose);

  # Compute internal actions
  ComputeInternalActions(S_obj, ext, iso_sol, dim, parse("verbose") = verbose);

  # Compute structure internal energy
  P := ComputePotentialEnergy(objs, parse("shear_contribution") = shear_contribution);

  # Solve hyperstatic equations
  hyper_eq := [];
  for i from 1 to nops(hyper_vars) do
    hyper_eq := [op(hyper_eq), diff(P, hyper_vars[i]) = hyper_disp[i]]
  end do;
  hyper_sol := op(solve(hyper_eq, hyper_vars));
  printf("DONE\n");

  return [op(hyper_sol), op(subs(hyper_sol, iso_sol))];
end proc: # HyperstaticSolver

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputePotentialEnergy := proc(
  objs::{ # Structure objects
    list({BEAM, ROD, SUPPORT, JOINT}),
    set( {BEAM, ROD, SUPPORT, JOINT})
  },
  {
    shear_contribution := false # Add shear contribution to the potential energy
  }, 
  $)

  description "Compute the internal potential energy of the structure";

  local i, P;

  P := 0;
  for i from 1 to nops(objs) do
    if IsBeam(objs[i]) or IsRod(objs[i]) then
      # Normal action N contribution
      if (member(N, map(lhs, objs[i][internal_actions]))) and
          (subs(objs[i][internal_actions](x), N(x)) <> 0) then
        P := P + integrate(
          subs(objs[i][internal_actions](x),
            N(x)^2/(2*objs[i][material][elastic_modulus]*objs[i][area])
          ), x = 0..objs[i][length]);
      end if;
      if shear_contribution then
        # Shear action Ty contribution
        if (member(Ty, map(lhs, objs[i][internal_actions]))) and
            (subs(objs[i][internal_actions](x), Ty(x)) <> 0) then
          P := P + integrate(
            subs(objs[i][internal_actions](x),
              objs[i][shear_stiff_factor][1]*Ty(x)^2/(2*objs[i][material][shear_modulus]*objs[i][area])
            ), x = 0..objs[i][length]);
        end if;
        # Shear action Tz contribution
        if (member(Tz, map(lhs, objs[i][internal_actions]))) and
            (subs(objs[i][internal_actions](x), Tz(x)) <> 0) then
          P := P + integrate(
            subs(objs[i][internal_actions](x),
              objs[i][shear_stiff_factor][2]*Tz(x)^2/(2*objs[i][material][shear_modulus]*objs[i][area])
            ), x = 0..objs[i][length]);
        end if;
      end if;
      # Bending moment action Mx contribution
      if (member(Mx, map(lhs, objs[i][internal_actions]))) and
          (subs(objs[i][internal_actions](x), Mx(x)) <> 0) then
        P := P + integrate(
          subs(objs[i][internal_actions](x),
            Mx(x)^2/(2*objs[i][material][shear_modulus]*objs[i][inertias][1])
          ), x = 0..objs[i][length]);
          end if;
      # Bending moment action My contribution
      if (member(My, map(lhs, objs[i][internal_actions]))) and
          (subs(objs[i][internal_actions](x), My(x)) <> 0) then
        P := P + integrate(
          subs(objs[i][internal_actions](x),
            My(x)^2/(2*objs[i][material][elastic_modulus]*objs[i][inertias][3])
          ), x = 0..objs[i][length]);
      end if;
      # Bending moment action Mz contribution
      if (member(Mz, map(lhs, objs[i][internal_actions]))) and
          (subs(objs[i][internal_actions](x), Mz(x)) <> 0) then
        P := P + integrate(
          subs(objs[i][internal_actions](x),
            Mz(x)^2/(2*objs[i][material][elastic_modulus]*objs[i][inertias][2])
          ), x = 0..objs[i][length]);
      end if;
    elif IsSupport(obj[i]) and IsCompliantSupport(obj[i]) then
      # Support reaction Fx contribution
      if (subs(objs[i][support_reactions], Fx) <> 0) and
          (objs[i][stiffness][1] <> infinity) then
        P := P + subs(objs[i][support_reactions], Fx^2/(2*objs[i][stiffness][1]));
      end if;
      # Support reaction Fy contribution
      if (subs(objs[i][support_reactions], Fy) <> 0) and
          (objs[i][stiffness][2] <> infinity) then
        P := P + subs(objs[i][support_reactions], Fy^2/(2*objs[i][stiffness][2]));
      end if;
      # Support reaction Fz contribution
      if (subs(objs[i][support_reactions], Fz) <> 0) and
          (objs[i][stiffness][3] <> infinity) then
        P := P + subs(objs[i][support_reactions], Fz^2/(2*objs[i][stiffness][3]));
      end if;
      # Support reaction Mx contribution
      if (subs(objs[i][support_reactions], Mx) <> 0) and
          (objs[i][stiffness][4] <> infinity) then
        P := P + subs(objs[i][support_reactions], Mx^2/(2*objs[i][stiffness][4]));
      end if;
      # Support reaction My contribution
      if (subs(objs[i][support_reactions], My) <> 0) and
          (objs[i][stiffness][5] <> infinity) then
        P := P + subs(objs[i][support_reactions], My^2/(2*objs[i][stiffness][5]));
      end if;
      # Support reaction Mz contribution
      if (subs(objs[i][support_reactions], Mz) <> 0) and
          (objs[i][stiffness][6] <> infinity) then
        P := P + subs(objs[i][support_reactions], Mz^2/(2*objs[i][stiffness][6]));
      end if;
    end if;
  end do;

  return P;
end proc: # PotentialEnergy

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

IsostaticSolver := proc(
  objs::{list(STRUCTURAL), set(STRUCTURAL)}, # Strucural objects
  ext::{list(LOAD), set(LOAD)},              # Load objects
  vars::list,                                # Variables to solve
  {
    dim::string      := "3D", # Dimension ("2D" or "3D")
    verbose::boolean := false # Verbose mode
  },
  $)

  description "Solve the isostatic structure equilibrium equation system given "
    "the structure objects, the external actions and the variables to solve";

  local eq, i, j, active_ext, sol, A, B, rank_eq, vars_tmp;

  # Compute structure equations
  printf("Message (in IsostaticSolver) computing the equilibrium equation for "
    "the isostatic structure...");
  eq := [];
    for i from 1 to nops(objs) do
    active_ext := {};
    for j from 1 to nops(ext) do
      if (ext[j][target] = objs[i][name]) then
        active_ext := {op(active_ext), ext[j]};
            end if;
        end do;
    eq := [op(eq), op(NewtonEuler(active_ext, objs[i], 0, parse("dim") = dim))];
    # Add joints and supports constraint equations
    if IsSupport(objs[i]) or IsJoint(objs[i]) then
      eq := [op(eq), op(objs[i][constraint_loads])];
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
  printf("DONE\n");

  # Structure equations check
  A, B    := LinearAlgebra[GenerateMatrix](eq, vars_tmp);
    rank_eq := LinearAlgebra[Rank](A);

  if (verbose) then
    printf("Message (in IsostaticSolver) structure equilibrium equations:\n");
        print(<op(eq)>);
    printf("Message (in IsostaticSolver) structure unknown variables:\n");
    print(vars_tmp);
    end if;

  if (rank_eq <> nops(vars_tmp)) then
    error "inconsistent system of equation, got %1 independent equations and "
      "%2 variables, check structure supports and joints",
      rank_eq, nops(vars_tmp);
  end if;

  # Solve structure equations
  printf("Message (in IsostaticSolver) computing the structure reaction forces... ");
  sol := simplify(op(solve(eq, vars_tmp)));
  printf("DONE\n");

  return sol;
end proc: # IsostaticSolver

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeInternalActions := proc(
  objs::{ # Structure objects
    list({BEAM, ROD}),
    set( {BEAM, ROD})
  },
  ext::{list(LOAD), set(LOAD)}, # Load objects
  sol::{list, set},             # Structure solution
  {
    dim::{string}      := "3D", # Dimension ("2D" or "3D")
    verbose::{boolean} := false # Verbose mode
  },
  $)

  description "Programmatic computation of internal actions for structure"
    "objects with given external actions and structure solution";

  local i, j, active_ext, subs_ext;

  # Substitute structure solution into loads
  subs_ext := map2(subs, sol, map(op,ext));

  for i from 1 to nops(objs) do
    # Extract active loads
    active_ext := {};
    for j from 1 to nops(subs_ext) do
      if (subs_ext[j][target] = objs[i][name]) then
        active_ext := {op(active_ext), subs_ext[j]};
      end if;
    end do;
    # Compute internal actions
    InternalActions(objs[i], active_ext, dim);
  end do;
end proc: # ComputeInternalActions

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

InternalActions := proc(
  obj::{BEAM, ROD}, # Structure object
  ext::{ # External actions
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  },
  {
    dim::{string} := "3D" # Dimension ("2D" or "3D")
  },
  $)

  description "Programmatic computation of internal actions for structure "
    "object with given external actions, it returns the internal actions as "
    "function of the axial variable 'x'";

  local i, ia, N_sol, Ty_sol, Tz_sol, Mx_sol, My_sol, Mz_sol;

  # 2D case
  if (dim = "2D") then

  N_sol  := 0;
  Ty_sol := 0;
  Mz_sol := 0;

  # Compute internal actions for concentrated loads as effect overlay
  for i from 1 to nops(ext) do
    if IsForce(ext[i]) then
      N_sol  := `simplify/piecewise`(N_sol  - piecewise(x >= ext[i][coordinate] and x <= obj[length], ext[i][components][1]), x);
      Ty_sol := `simplify/piecewise`(Ty_sol + piecewise(x >= ext[i][coordinate] and x <= obj[length], ext[i][components][2]), x);
      Mz_sol := `simplify/piecewise`(Mz_sol + integrate(piecewise(x >= ext[i][coordinate] and x <= obj[length], ext[i][components][2]), x = 0..x), x);
    elif IsMoment(ext[i]) then
      Mz_sol := `simplify/piecewise`(Mz_sol - piecewise(x >= ext[i][coordinate] and x <= obj[length], ext[i][components][3]), x);
    elif IsQForce(ext[i]) then
      N_sol  := `simplify/piecewise`(N_sol  - integrate(ext[i][components][1](x), x = 0..x), x);
      Ty_sol := `simplify/piecewise`(Ty_sol + integrate(ext[i][components][2](x), x = 0..x), x);
      Mz_sol := `simplify/piecewise`(Mz_sol + integrate(integrate(ext[i][components][2](x), x = 0..x), x = 0..x), x);
    elif IsQMoment(FMQ[i]) then
      Mz_sol := `simplify/piecewise`(Mz_sol - integrate(ext[i][components][3](x), x = 0..x), x);
    end if;
    end do;

  ia := [
    N = unapply(N_sol, x), Ty = unapply(Ty_sol, x), Mz = unapply(Mz_sol, x)
            ];

  # 3D case
  elif (dim = "3D") then

    N_sol  := 0;
    Ty_sol := 0;
    Tz_sol := 0;
    Mx_sol := 0;
    My_sol := 0;
    Mz_sol := 0;

    # Compute internal actions for concentrated loads as effect overlay
    for i from 1 to nops(ext) do
      if IsForce(ext[i]) then
        N_sol  := `simplify/piecewise`(N_sol  - piecewise(x >= ext[i][coordinate] and x <= obj[length], ext[i][components][1]), x);
        Ty_sol := `simplify/piecewise`(Ty_sol + piecewise(x >= ext[i][coordinate] and x <= obj[length], ext[i][components][2]), x);
        Tz_sol := `simplify/piecewise`(Tz_sol + piecewise(x >= ext[i][coordinate] and x <= obj[length], ext[i][components][3]), x);
        My_sol := `simplify/piecewise`(My_sol + integrate(piecewise(x >= ext[i][coordinate] and x <= obj[length], ext[i][components][3]), x = 0..x), x);
        Mz_sol := `simplify/piecewise`(Mz_sol + integrate(piecewise(x >= ext[i][coordinate] and x <= obj[length], ext[i][components][2]), x = 0..x), x);
      elif IsMoment(ext[i]) then
        Mx_sol := `simplify/piecewise`(Mx_sol - piecewise(x >= ext[i][coordinate] and x <= obj[length], ext[i][components][1]), x);
        My_sol := `simplify/piecewise`(My_sol + piecewise(x >= ext[i][coordinate] and x <= obj[length], ext[i][components][2]), x);
        Mz_sol := `simplify/piecewise`(Mz_sol - piecewise(x >= ext[i][coordinate] and x <= obj[length], ext[i][components][3]), x);
      elif IsQForce(ext[i]) then
        N_sol  := `simplify/piecewise`(N_sol  - integrate(ext[i][components][1](x), x = 0..x), x);
        Ty_sol := `simplify/piecewise`(Ty_sol + integrate(ext[i][components][2](x), x = 0..x), x);
        Tz_sol := `simplify/piecewise`(Tz_sol + integrate(ext[i][components][3](x), x = 0..x), x);
        My_sol := `simplify/piecewise`(My_sol + integrate(integrate(ext[i][components][3](x), x = 0..x), x = 0..x), x);
        Mz_sol := `simplify/piecewise`(Mz_sol + integrate(integrate(ext[i][components][2](x), x = 0..x), x = 0..x), x);
      elif IsQMoment(FMQ[i]) then
        Mx_sol := `simplify/piecewise`(Mx_sol - integrate(ext[i][components][1](x), x = 0..x), x);
        My_sol := `simplify/piecewise`(My_sol + integrate(ext[i][components][2](x), x = 0..x), x);
        Mz_sol := `simplify/piecewise`(Mz_sol - integrate(ext[i][components][3](x), x = 0..x), x);
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
    "Message (in InternalActions) updating %s %s's internal actions... ",
    obj[type], obj[name]
    );
  obj[internal_actions] := ia;
  printf("DONE\n");

  return ``;
end proc: # InternalActions

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ComputeDisplacements := proc(
  objs::{ # Structure objects
    list({BEAM, ROD, SUPPORT, JOINT}),
    set( {BEAM, ROD, SUPPORT, JOINT})
  },
  ext::{ # External actions
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  },
  sol::{list, set}, # Solution of the structure
  {
    dim::{string} := "3D" # Dimension ("2D" or "3D")
    shear_contribution := false # Add shear contribution to the potential energy
  },
  $)

  description "Compute the Structure displacements";

  local dummy_Fx, dummy_Fy, dummy_Fz, dummy_Mx, dummy_My, dummy_Mz, i, obj_copy, dummy_loads, subs_null_dummy;

  # Cicle on the structure objects
  for i from 1 to nops(objs) do
    # Create a copy of the object
    obj_copy := copy(objs[i]);

    # Beam
    if IsBeam(obj_copy) then
      # Create dummy loads 
      dummy_Fx := MakeForce(['dFx',0,0], 'x', {obj_copy}, obj_copy[frame]);
      dummy_Fy := MakeForce([0,'dFy',0], 'x', {obj_copy}, obj_copy[frame]);
      dummy_Fz := MakeForce([0,0,'dFz'], 'x', {obj_copy}, obj_copy[frame]);
      dummy_Mx := MakeMoment(['dMx',0,0], 'x', {obj_copy}, obj_copy[frame]);
      dummy_My := MakeMoment([0,'dMy',0], 'x', {obj_copy}, obj_copy[frame]);
      dummy_Mz := MakeMoment([0,0,'dMz'], 'x', {obj_copy}, obj_copy[frame]);

      dummy_loads := {dummy_Fx, dummy_Fy, dummy_Fz, dummy_Mx, dummy_My, dummy_Mz};

      # Compute internal actions of the object copy
      ComputeInternalActions(obj_copy, ext union dummy_loads, sol, dim, parse("verbose") = verbose);

      # Compute object potential energy
      P := ComputePotentialEnergy(obj_copy, parse("shear_contribution") = shear_contribution);

      # null dummy loads substitution list
      subs_null_dummy := ['dFx', 'dFy', 'dFz', 'dMx', 'dMy', 'dMz'] =~ [0,0,0,0,0,0];

      # Compute displacements
      objs[i][displacement][1] := subs(subs_null_dummy, diff(P, dFx));
      objs[i][displacement][2] := subs(subs_null_dummy, diff(P, dFy));
      objs[i][displacement][3] := subs(subs_null_dummy, diff(P, dFz));
      objs[i][displacement][4] := subs(subs_null_dummy, diff(P, dMx));
      objs[i][displacement][5] := subs(subs_null_dummy, diff(P, dMy));
      objs[i][displacement][6] := subs(subs_null_dummy, diff(P, dMz));

    # Rod
    elif IsRod(obj_copy) then
      # Create dummy loads
      dummy_Fx := MakeForce(['dFx',0,0], 'x', {obj_copy}, obj_copy[frame]);

      dummy_loads := {dummy_Fx};

      # Compute internal actions of the object copy
      ComputeInternalActions(obj_copy, ext union dummy_loads, sol, dim, parse("verbose") = verbose);

      # Compute object potential energy
      P := ComputePotentialEnergy(obj_copy, parse("shear_contribution") = shear_contribution);

      # null dummy loads substitution list
      subs_null_dummy := ['dFx'] =~ [0];

      # Compute displacements
      objs[i][displacement][1] := subs(subs_null_dummy, diff(P, dFx));

    # Support
    elif IsSupport(obj_copy) then 
      # Create dummy loads on constrained directions and add to support reactions
      if (obj_copy[constrained_dof][1] <> 0) then
        obj_copy[support_reactions][1] := obj_copy[support_reactions][1] + 'dFx'; 
      end if;
      if (obj_copy[constrained_dof][2] <> 0) then
        obj_copy[support_reactions][2] := obj_copy[support_reactions][2] + 'dFy'; 
      end if;
      if (obj_copy[constrained_dof][3] <> 0) then
        obj_copy[support_reactions][3] := obj_copy[support_reactions][3] + 'dFz'; 
      end if;
      if (obj_copy[constrained_dof][4] <> 0) then
        obj_copy[support_reactions][4] := obj_copy[support_reactions][4] + 'dMx'; 
      end if;
      if (obj_copy[constrained_dof][5] <> 0) then
        obj_copy[support_reactions][5] := obj_copy[support_reactions][5] + 'dMy'; 
      end if;
      if (obj_copy[constrained_dof][6] <> 0) then
        obj_copy[support_reactions][6] := obj_copy[support_reactions][6] + 'dMz'; 
      end if;

      # Compute potential energy of the object copy
      P := ComputePotentialEnergy(obj_copy);

      # null dummy loads substitution list
      subs_null_dummy := ['dFx', 'dFy', 'dFz', 'dMx', 'dMy', 'dMz'] =~ [0,0,0,0,0,0];

      # Compute displacements
      objs[i][displacement][1] := subs(subs_null_dummy, diff(P, dFx));
      objs[i][displacement][2] := subs(subs_null_dummy, diff(P, dFy));
      objs[i][displacement][3] := subs(subs_null_dummy, diff(P, dFz));
      objs[i][displacement][4] := subs(subs_null_dummy, diff(P, dMx));
      objs[i][displacement][5] := subs(subs_null_dummy, diff(P, dMy));
      objs[i][displacement][6] := subs(subs_null_dummy, diff(P, dMz));
    end if;
  end do;

  return ``;
end proc; # ComputeDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module:

# That's all folks!
