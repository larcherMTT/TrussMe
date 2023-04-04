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

  # TODO: Hide module content with dummy procedures as:
  #   procname := proc(inputs); _procname(inputs); end proc;
  # so the source code is not visible by showstat and showsource commands.

  global  ground;

  local m_LAST;
  local m_LEM;
  local m_earth;
  local m_gravity;
  local m_VerboseMode;
  local m_WarningMode;
  local m_TimeLimit;
  local m_KeepVeiled;
  local m_BeamColor;
  local m_RodColor;
  local m_RigidBodyColor;
  local m_CompliantSupportColor;
  local m_SupportColor;
  local m_CompliantJointColor;
  local m_JointColor;
  local m_EarthColor;
  local m_StoredData;

  option  package,
          load   = ModuleLoad,
          unload = ModuleUnload;

  description "A Maple Library for Truss Elements Structures.";

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Info := proc() # REVIEWED

    description "Print 'TrussMe' module information.";

    printf(
      "+--------------------------------------------------------------------------+\n"
      "| 'TrussMe' module version 0.0 - BSD 3-Clause License - Copyright (c) 2023  |\n"
      "| Current version authors:                                                  |\n"
      "|   Matteo Larcher and Davide Stocco.                                       |\n"
      "+--------------------------------------------------------------------------+\n"
    );
    return NULL;
  end proc: # Info

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleLoad := proc() # REVIEWED

    description "Module 'TrussMe' load procedure.";

    local lib_base_path, i;

    lib_base_path := NULL;
    for i in [libname] do
      if (StringTools:-Search("TrussMe", i) <> 0) then
        lib_base_path := i;
      end if;
    end do;
    if (lib_base_path = NULL) then
      error "cannot find 'TrussMe' library.";
    end if;

    TrussMe:-TypeRegister();
    TrussMe:-InitTrussMe();
    TrussMe:-Protect();
    return NULL;
  end proc: # ModuleLoad

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleUnload := proc() # REVIEWED

    description "Module 'TrussMe' unload procedure.";

    m_VerboseMode := 1;
    m_WarningMode := true;
    TrussMe:-UnAssignData();
    TrussMe:-Unprotect();
    return NULL;
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export TypeRegister := proc() # REVIEWED

    description "Register 'TrussMe' module types.";

    TypeTools:-AddType(EARTH, IsEarth);
    TypeTools:-AddType(FRAME, IsFrame);
    TypeTools:-AddType(POINT, IsPoint);
    TypeTools:-AddType(VECTOR, IsVector);
    TypeTools:-AddType(BEAM, IsBeam);
    TypeTools:-AddType(ROD, IsRod);
    TypeTools:-AddType(RIGID_BODY, IsRigidBody);
    TypeTools:-AddType(FORCE, IsForce);
    TypeTools:-AddType(MOMENT, IsMoment);
    TypeTools:-AddType(QFORCE, IsQForce);
    TypeTools:-AddType(QMOMENT, IsQMoment);
    TypeTools:-AddType(SUPPORT, IsSupport);
    TypeTools:-AddType(JOINT, IsJoint);
    TypeTools:-AddType(MATERIAL, IsMaterial);
    TypeTools:-AddType(STRUCTURE, IsStructure);
    return NULL;
  end proc: # TypeRegister

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export InitTrussMe := proc() # REVIEWED

    description "Initialize 'TrussMe' module internal variables.";

    TrussMe:-InitLAST();

    ground := <<1, 0, 0, 0>|
                        <0, 1, 0, 0>|
                        <0, 0, 1, 0>|
                        <0, 0, 0, 1>>;

    m_gravity := [0, 0, 0];

    m_VerboseMode           := 1;
    m_WarningMode           := true;
    m_TimeLimit             := 5;
    m_BeamColor             := "SteelBlue";
    m_RodColor              := "Niagara DarkOrchid";
    m_RigidBodyColor        := "Indigo";
    m_CompliantSupportColor := "DarkGreen";
    m_SupportColor          := "DarkOrange";
    m_CompliantJointColor   := "LightSalmon";
    m_JointColor            := "MediumSeaGreen";
    m_EarthColor            := "Firebrick";
    m_StoredData            := [];

    m_earth := table({
      "type"             = EARTH,
      "name"             = "earth",
      "length"           = 0,
      "frame"            = ground,
      "admissible_loads" = [1, 1, 1, 1, 1, 1]
      }):

    return NULL;
  end proc: # InitTrussMe

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Protect := proc() # REVIEWED

    description "Protect 'TrussMe' module variables.";

    protect(
      # Global variables
      'ground',
      # Types
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
    return NULL;
  end proc: # Protect

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Unprotect := proc() # REVIEWED

    description "Unprotect 'TrussMe' module variables.";

    unprotect(
      # Global variables
      'ground',
      # Types
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
    return NULL;
  end proc: # Unprotect

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export CheckInit := proc( $ ) # REVIEWED

    description "Check if the 'LAST' object is initialized.";

    if (m_LAST = NULL) then
      error "the 'LAST' object is not initialized, use 'TrussMe:-InitLAST(...)' "
        "or other appropriate initialization methods first.";
    end if;
    if (m_LEM = NULL) then
      error "the 'LEM' object is not initialized, use 'TrussMe:-InitLAST(...)' "
        "or other appropriate initialization methods first.";
    end if;
    return NULL;
  end proc: # CheckInit

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export InitLAST := proc( # REVIEWED
    label::{symbol, string} := NULL,
    $)

    description "Initialize the 'LAST' object with veiling label <label>.";

    m_LAST := Object(LAST);
    m_LAST:-InitLEM(m_LAST, label);
    m_LEM  := m_LAST:-GetLEM(m_LAST);
    return NULL;
  end proc: # InitLAST

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ClearLAST := proc( $ )

    description "Clear the 'LAST' (and 'LEM') object.";

    m_LAST := NULL;
    m_LEM  := NULL;
    return NULL;
  end proc: # ClearLAST

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetLAST := proc( # REVIEWED
    obj::LAST,
    $)

    description "Set the 'LAST' (and 'LEM') object <obj>.";

    m_LAST := obj;
    m_LEM  := m_LAST:-GetLEM(m_LAST);
    return NULL;
  end proc: # SetLAST

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetLAST := proc( $ )::LAST; # REVIEWED

    description "Get the 'LAST' object.";

    return m_LAST;
  end proc: # GetLAST

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetLEM := proc( # REVIEWED
    obj::LEM,
    $)

    description "Set the 'LEM' object <obj>.";

    m_LAST:-SetLEM(m_LAST, obj);
    m_LEM := obj;
    return NULL;
  end proc: # SetLEM

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetLEM::static := proc( $ )::LEM;

    description "Get the 'LEM' object.";

    return m_LEM;
  end proc: # GetLEM

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export `union` := proc( # REVIEWED
  A::{list, set},
  B::{list, set},
  $)::{list, set};

  option overload;

  description "Extension of union operator to list or set objects <A> and <B>.";

  if type(A, 'set') and type(B, 'set') then
    return {op(A), op(B)};
  else
    return [op(A), op(B)];
  end if;
end proc: # union

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export SetModuleOptions := proc( # REVIEWED
  {
    VerboseMode::{integer, nothing} := NULL,
    WarningMode::{boolean, nothing} := NULL,
    TimeLimit::{constant, nothing}  := NULL
  },
  $)

  description "Set the module options: <VerboseMode> = [0, 1, 2], <WarningMode> "
    "= [true, false] and <TimeLimit> = [0, +inf].";

  if (VerboseMode <> NULL) then
    if (VerboseMode < 0) or (VerboseMode > 2) then
      error "invalid verbose mode detected.";
    else
      m_VerboseMode := VerboseMode;
    end if;
  end if;

  if (WarningMode <> NULL) then
    if not type(WarningMode, boolean) then
      error "invalid warning mode detected.";
    else
      m_WarningMode := WarningMode;
    end if;
  end if;

  if (TimeLimit <> NULL) then
    if (TimeLimit < 0) then
      error "invalid time limit detected.";
    else
      m_TimeLimit := TimeLimit;
    end if;
  end if;
  return NULL;
end proc: # SetModuleOptions

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export EnableVerboseMode::static := proc( $ ) # REVIEWED

  description "Enable the verbose mode of the module.";

  m_VerboseMode := 1;
  return NULL;
end proc: # EnableVerboseMode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DisableVerboseMode::static := proc( $ ) # REVIEWED

  description "Disable the verbose mode of the module.";

  m_VerboseMode := false;
  return NULL;
end proc: # DisableVerboseMode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export EnableWarningMode::static := proc( $ ) # REVIEWED

  description "Enable the warning mode of the module.";

  m_WarningMode := true;
  return NULL;
end proc: # EnableWarningMode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DisableWarningMode::static := proc( $ ) # REVIEWED

  description "Disable the warning mode of the module.";

  m_WarningMode := 0;
  return NULL;
end proc: # DisableWarningMode

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsEarth := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Test if an object <obj> is the EARTH object.";

  return evalb(obj["type"] = EARTH);
end proc: # IsEarth

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export SetGravity := proc( # REVIEWED
  vec::{list, Vector},
  $)

  description "Set gravity vector with [x, y, z]^T components of <vec>.";

  if type(vec, list) and (nops(vec) = 3) then
    m_gravity := <op(vec)>;
  elif type(vec, Vector) and (nops(vec) = 3) then
    m_gravity := vec;
  else
    error "invalid gravity vector detected.";
  end if;
  return NULL;
end proc: # SetGravity

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetGravity := proc( $ )::Vector; # REVIEWED

  description "Get gravity vector.";

  return m_gravity;
end proc: # GetGravity

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Norm2 := proc( # REVIEWED
  vec::{list, vector},
  $)::algebraic;

  description "Compute the Euclidean norm of the input vector <vec>.";

  return sqrt(add(x, x in vec^~2));
end proc: # Norm2

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ListPadding := proc( # REVIEWED
  lst::{algebraic, list, Vector},
  n::integer,
  value::algebraic := 0,
  $)::{list, Vector};

  description "Pad a list or vector <lst> with <value> to have <n> elements.";

  local out;

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
  return out;
end proc: # ListPadding

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Show := proc( # REVIEWED
  tab::table,
  $)

  description "Show the content of a table <tab>.";

  print(tab = tab["type"](op(op(tab))));
  return NULL;
end proc: # Show

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetNames := proc( # REVIEWED
  objs::{
    list({MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH}),
    set({MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH})
  }, $)::{list({string}), set({string})};

  description "Get names of a list or set of structural objects <objs>.";

  if type(objs, 'set') then
    return {seq(objs[i]["name"], i = 1..nops(objs))};
  else
    return [seq(objs[i]["name"], i = 1..nops(objs))];
  end if;
end proc: # GetNames

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetObjByName := proc( # REVIEWED
  name::string,
  objs::{
    list({MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH}),
    set({MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH})
  }, $)::{anything, MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH};

  description "Get object which name field is <name> from a list or set of "
    "objects <objs>.";

  local out, obj;

  out := NULL;
  for obj in objs do
    if (obj["name"] = name) then
      out := eval(obj); # Do not remove eval
      break;
    end if;
  end do;
  return eval(out);
end proc: # GetObjByName

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export GetObjsByType := proc( # REVIEWED
  types::{list(symbol), set(symbol)},
  objs::{list, set},
  $)::list;

  description "Get objects which type field is in <types> from a list or set of "
    "objects <objs>.";

  local out, obj;

  out := [];
  for obj in objs do
    if (obj::convert(types, set)) then
      out := out union [eval(obj)]; # Do not remove eval
    end if;
  end do;
  return eval(out);
end proc: # GetObjsByType

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Simplify := proc( # REVIEWED
  obj::anything,
  opt::anything := NULL,
  $)::anything;

  description "Simplify an algebraic expression <obj>.";

  local out, time_limit;

  time_limit := `if`(procname::indexed, op(procname), m_TimeLimit); # FIXME unclear
  try
    out := timelimit(time_limit, simplify(obj, opt));
  catch:
    WARNING("time limit of %1s exceeded for simplify operation, raw solutions "
      "is returned. The input <time_limit> can be modified by setting it in the "
      "TrussMe:-SetModuleOptions(...) procedure.", time_limit
    );
    out := obj;
  end try:
  return out;
end proc: # Simplify

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export AssignData := proc( # REVIEWED
  data::{list, set},
  $)

  description "Assign the list <x> to the local variable <m_StoredData>.";

  m_StoredData := data;
  return NULL;
end proc: # AssignData

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export UnAssignData := proc( # REVIEWED
  $)

  description "Unassign the local variable <m_StoredData>.";

  m_StoredData := [];
  return NULL;
end proc: # UnAssignData

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Subs := proc( # REVIEWED
  # _passed
  )::anything;

  description "Perform subs command neglecting sub-lists and sub-sets from the "
    "substitution list.";

  map(x -> map(remove, y -> type(y, {list, set}), x), [_passed[1..-2]]);
  return subs(op(%), _passed[-1]);
end proc; # Subs

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Diff := proc( # REVIEWED
  # _passed
  {
  veils := NULL
  })::anything;

  description "Perform diff command on veiled expressions given the veiling "
    "list <veils>.";

  # TODO: optimize this procedure

  local subs_diff, d_vars, v2f, f2v, last, veils_copy;

  if (veils = NULL) then
    veils_copy := [];
    last := -1
  else
    veils_copy := veils;
    last := -2;
  end if;

  # Get the variables to be differentiated
  d_vars := _passed[2..last];

  # Veil to functions substitution list
  v2f := map(x -> lhs(x) =~ lhs(x)(d_vars), veils_copy);

  # Function to veils substitution list
  f2v := rhs~(v2f) =~ lhs~(v2f);

  subs(v2f, veils_copy);
  diff(lhs~(%), d_vars) =~ TrussMe:-Simplify(diff(rhs~(%), d_vars));
  subs_diff := lhs~(%) =~ TrussMe:-Simplify(subs(op(ListTools:-Reverse(%)),rhs~(%))):

  # Compute the derivative of the veiled expression
  subs(subs_diff, diff(subs(v2f, _passed[1]), d_vars));

  # Substitute back the veils
  return subs(f2v, %);
end proc; # Diff

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export InverseFrame := proc( # REVIEWED
  RF::FRAME,
  $)::FRAME;

  description "Inverse transformation matrix of an affine transformation <RF>.";

  LinearAlgebra:-Transpose(RF[1..3, 1..3]);
  return <<% | -% . RF[1..3, 4]>,
          <0 | 0 | 0 | 1>>;
end proc: # InverseFrame

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsFrame := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Check if the input object <obj> is a FRAME object.";

  if (type(obj, Matrix)) and
     (LinearAlgebra:-RowDimension(obj) = 4) and
     (LinearAlgebra:-ColumnDimension(obj) = 4) then
    return true;
  else
    return false;
  end if;
end proc: # IsFrame

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Rotate := proc( # REVIEWED
  axis::{symbol, string},
  angle::algebraic,
  $)::FRAME;

  description "Transformation matrix corresponding to the rotation <angle> "
    "around the given <axis>";

  if (axis = 'X') or (axis = "X") then
    return <<1, 0,           0,          0>|
            <0, cos(angle),  sin(angle), 0>|
            <0, -sin(angle), cos(angle), 0>|
            <0, 0,           0,          1>>;
  elif (axis = 'Y') or (axis = "Y") then
    return <<cos(angle), 0, -sin(angle), 0>|
            <0,          1, 0,           0>|
            <sin(angle), 0, cos(angle),  0>|
            <0,          0, 0,           1>>;
  elif (axis = 'Z') or (axis = "Z") then
    return <<cos(angle),  sin(angle), 0, 0>|
            <-sin(angle), cos(angle), 0, 0>|
            <0,           0,          1, 0>|
            <0,           0,          0, 1>>;
  else
    error "invalid axis detected.";
  end if;
end proc: # Rotate

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Translate := proc( # REVIEWED
  x::algebraic,
  y::algebraic,
  z::algebraic,
  $)::FRAME;

  description "Transformation matrix corresponding to the translation <x, y, z>.";

  return <<1, 0, 0, 0>|
          <0, 1, 0, 0>|
          <0, 0, 1, 0>|
          <x, y, z, 1>>;
end proc: # Translate

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Origin := proc( # REVIEWED
  RF::FRAME,
  $)::Vector;

  description "Extract the origin of the reference frame <RF>.";

  return <RF[1, 4], RF[2, 4], RF[3, 4], 1>;
end proc: # Origin

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsPoint := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Check if the input object <obj> is a POINT object.";

  if (type(obj, list)) and
     (nops(obj) = 3) then
    return true;
  else
    return false;
  end if;
end proc: # IsPoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsVector := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Check if the input object <obj> is a VECTOR object.";

  if (type(obj, list)) and
     (nops(obj) = 3) then
    return true;
  else
    return false;
  end if;
end proc: # IsVector

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Uvec := proc( # REVIEWED
  axis::symbol,
  RF::FRAME := ground,
  $)::Vector;

  description "Extract the unit vector of the reference frame <RF> along the "
    "given <axis>.";

  if (axis = 'X') then
    return <RF[1, 1], RF[2, 1], RF[3, 1], 0>;
  elif (axis = 'Y') then
    return <RF[1, 2], RF[2, 2], RF[3, 2], 0>;
  elif (axis = 'Z') then
    return <RF[1, 3], RF[2, 3], RF[3, 3], 0>;
  else
    error "invalid axis detected";
  end if;
end proc: # Uvec

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export UvecX := proc( # REVIEWED
  RF::FRAME := ground,
  $)::Vector;

  description "Extract the x-axis unit vector of the reference frame <RF>.";

  return <RF[1, 1], RF[2, 1], RF[3, 1], 0>;
end proc: # UvecX

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export UvecY := proc( # REVIEWED
  RF::FRAME := ground,
  $)::Vector;

  description "Extract the y-axis unit vector of the reference frame <RF>.";

  return <RF[1, 2], RF[2, 2], RF[3, 2], 0>;
end proc: # UvecY

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export UvecZ := proc( # REVIEWED
  RF::FRAME := ground,
  $)::Vector;

  description "Extract the z-axis unit vector of the reference frame <RF>.";

  return <RF[1, 3], RF[2, 3], RF[3, 3], 0>;
end proc: # UvecZ

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export Project := proc( # REVIEWED
  x::{list, Vector},
  RF_ini::FRAME,
  RF_end::FRAME,
  $)::{list, Vector};

  description "Project <x,y,z>, or vector <x,y,z,0>, or point <x,y,z,1> from "
    "reference frame <RF_ini> to reference frame <RF_end>.";

  local x_tmp, out, i;

  # Pad input vector with 0 if its length is 3
  if not (nops(x) = 3) and
     not (nops(x) = 4) then
    error "invalid input vector/point <x> detected";
  end if;
  x_tmp := TrussMe:-ListPadding(convert(x, Vector), 4);

  # Try to compare RF_end and RF_ini
  try
    # FIXME: problems with floats (floats not handled error)
    map(evalb, evala(simplify(RF_end)) =~ evala(simplify(RF_ini)));
  catch:
    map(evalb, RF_end =~ RF_ini);
  end try;

  if has(%, false) then
    TrussMe:-InverseFrame(RF_end).RF_ini.x_tmp;
    out := TrussMe:-Simplify([seq(%[i], i = 1..nops(x))]);
  else
    out := x_tmp[1..nops(x)];
  end if;

  return convert(out, whattype(x));
end proc: # Project

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeMaterial := proc({ # REVIEWED
    material_name::string      := "DeafultSteel",
    elastic_modulus::algebraic := 210.0E+09,
    poisson_ratio::algebraic   := 0.3,
    shear_modulus::algebraic   := elastic_modulus/(2*(1+poisson_ratio)),
    density::algebraic         := 7.4E+03
  }, $)::MATERIAL;

  description "Define a MATERIAL object with inputs: name of the material "
    "<material_name>, elastic modulus <elastic_modulus> (default = 210.0E9 Pa), "
    "Poisson ratio <poisson_ratio> (default = 0.3), shear modulus <shear_modulus> "
    "(default = E/(2*(1+nu))), density <density> (default = 7.4E3 kg/m^3).";

  return table({
    "type"            = MATERIAL,
    "name"            = material_name,
    "elastic_modulus" = elastic_modulus,
    "poisson_ratio"   = poisson_ratio,
    "shear_modulus"   = shear_modulus,
    "density"         = density
    });
end proc: # DefineMaterial

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsMaterial := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Check if the input object <obj> is a MATERIAL object.";

  return evalb(obj["type"] = MATERIAL);
end proc: # IsMaterial

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeForce := proc( # REVIEWED
  components::{list, Vector},
  coords::{algebraic, list(algebraic)},
  obj::{BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH},
  RF::FRAME := ground,
  $)::FORCE;

  description "Define a 'FORCE' object with inputs: force components <components>, "
    "force application axial coordinate <coords>, target object <obj>, and "
    "optional reference frame <RF> in which the force is defined (default = "
    "ground).";

  local proj_components, admissible_components;

  # Check input argument
  if not (nops(components) = 3) then
    error "invalid input vector <components> detected.";
  end if;

  proj_components       := TrussMe:-Project(components, RF, obj["frame"]);
  admissible_components := convert(proj_components .~ <obj["admissible_loads"][1..3]>, list);
  if (proj_components <> admissible_components) and m_WarningMode then
    ["x_comp", "y_comp", "z_comp"] =~ convert(proj_components .~ <eval(
        map((x -> evalb(x = 0)), obj["admissible_loads"][1..3]), [true = 1, false = 0]
      )>, list);
    WARNING("Force components are not admissible for the target object. The "
      "following components will be ignored: %1", remove(x -> rhs(x) = 0, %));
  end if;

  if TrussMe:-IsSupport(obj) or TrussMe:-IsJoint(obj) then
    if (TrussMe:-ListPadding(coords, 3) <> [0, 0, 0]) then
      error "only null axial coordinate is accepted for 'SUPPORT' and 'JOINT' "
        "type objects";
    end if;
  end if;

  return table({
    "type"       = FORCE,
    "components" = admissible_components,
    "coordinate" = TrussMe:-ListPadding(coords, 3),
    "target"     = obj["name"]
    });
end proc: # MakeForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsForce := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Check if the object <obj> is a FORCE object";

  return evalb(obj["type"] = FORCE);
end proc: # IsForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeMoment := proc( # REVIEWED
  components::{list, Vector},
  coords::{algebraic, list(algebraic)},
  obj::{BEAM, RIGID_BODY, SUPPORT, JOINT, EARTH},
  RF::FRAME := ground,
  $)::MOMENT;

  description "Define a 'MOMENT' object with inputs: moment components <components>, "
    "moment application axial coordinate <coords>, target object <obj>, and "
    "optional reference frame <RF> in which the moment is  defined (default = "
    "ground).";

  local proj_components, admissible_components, out;

  # Check input argument
  if not (nops(components) = 3) then
    error "invalid input vector <components> detected.";
  end if;

  proj_components       := TrussMe:-Project(components, RF, obj["frame"]);
  admissible_components := convert(proj_components .~ <obj["admissible_loads"][4..6]>, list);
  if (proj_components <> admissible_components) and m_WarningMode then
    ["x_comp", "y_comp", "z_comp"] =~ convert(proj_components .~ <eval(
        map((x -> evalb(x = 0)), obj["admissible_loads"][4..6]), [true = 1, false = 0]
      )>, list);
    WARNING("Moment components are not admissible for the target object. The "
      "following components will be ignored: %1", remove(x -> rhs(x) = 0, %));
  end if;

  if (TrussMe:-IsSupport(obj) or TrussMe:-IsJoint(obj)) and
     (TrussMe:-ListPadding(coords, 3) <> [0, 0, 0]) then
    error "only null axial coordinate is accepted for 'SUPPORT' and 'JOINT' "
      "type objects";
  end if;

  return table({
    "type"       = MOMENT,
    "components" = admissible_components,
    "coordinate" = TrussMe:-ListPadding(coords, 3),
    "target"     = obj["name"]
    });
end proc: # MakeMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsMoment := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Check if the object <obj> is a MOMENT object.";

  return evalb(obj["type"] = MOMENT);
end proc: # IsMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeQForce := proc( # REVIEWED
  components::{procedure, list(algebraic)},
  obj::{BEAM, ROD},
  RF::FRAME := ground,
  {
    ell_min::algebraic := 0,
    ell_max::algebraic := obj["length"]
  }, $)::QFORCE;

  description "Define a QFORCE object with inputs: distributed load target "
    "object components <components>, target object <obj>, optional reference "
    "frame <RF> in which the load components are defined (default = ground), "
    "and optional initial <ell_min> and final <ell_max> application points "
    "(axial coordinates).";

  local proj_components, x;

  if type(components, procedure) then
    proj_components := unapply(
      TrussMe:-Project(components(x), RF, obj["frame"]), x
    );
  else
    proj_components := x -> piecewise(
      (ell_min <= x) and (x <= ell_max), TrussMe:-Project(components, RF, obj["frame"]), 0
    );
  end if;

  if TrussMe:-IsRod(obj) then
    if (proj_components(x)[2] <> 0) or (proj_components(x)[3] <> 0) then
      error "only axial loads are accepted in ROD objects"
    end if;
  end if;

  return table({
    "type"        = QFORCE,
    "components"  = proj_components,
    "target"      = obj["name"]
    });
end proc: # MakeQForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsQForce := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Check if the object <obj> is a QFORCE object.";

  return evalb(obj["type"] = QFORCE);
end proc: # IsQForce

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeQMoment := proc(
  components::{procedure, list(algebraic)},
  obj::BEAM,
  RF::FRAME := ground,
  {
    ell_min::algebraic := 0,
    ell_max::algebraic := obj["length"]
  }, $)::QMOMENT;

  description "Define a QMOMENT object with inputs: distributed torque target "
    "object components <components>, target object <obj>, optional reference "
    "frame <RF> in which the load components are defined (default = ground), "
    "and optional initial <ell_min> and final <ell_max> application points "
    "(axial coordinates).";

  local proj_components, x;

  if type(components, procedure) then
    proj_components := unapply(TrussMe:-Project(components(x), RF, obj["frame"]), x);
  else
    proj_components := x -> piecewise(
      (ell_min <= x) and (x <= ell_max), TrussMe:-Project(components, RF, obj["frame"]), 0
    );
  end if;

  return table({
    "type"        = QMOMENT,
    "components"  = proj_components,
    "target"      = obj["name"]
    });
end proc: # MakeQMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsQMoment := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Check if the object <obj> is a QMOMENT object";

  return evalb(obj["type"] = QMOMENT);
end proc: # IsQMoment

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeSupport := proc(
  name::string,                                        # Support name
  constrained_dof::list,                               # Constrained degree of freedom
  objs::{list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})}, # Target objects
  coords::list,                                        # Support locations
  RF::FRAME := ground,                                 # Reference frame of the support
  {
    stiffness::{procedure,list(algebraic)} := [ # Stiffness components (default = infinite)
      infinity, infinity, infinity,
      infinity, infinity, infinity
    ] *~ constrained_dof
  }, $)::SUPPORT;

  description "Make a SUPPORT object with inputs: support name <name>, constrained "
    "degrees of freedom <constrained_dof>, target objects <objs>, support locations "
    "<coords>, and optional reference frame <RF> in which the support is defined "
    "(default = ground). The optional input <stiffness> is a list of stiffness "
    "components (default = infinite) in the order: [ktx, kty, ktz, krx, kry, krz].";

  local S, J_tmp, i, j, sr_F_names, sr_F_values_tmp, sr_M_names, sr_M_values_tmp,
    S_stiffness, x, obj_coords;

  # Substitute -1 entries of coords with the corresponding object length
  obj_coords := [seq(`if`(coords[i] = -1, objs[i]["length"], coords[i]), i = 1..nops(coords))];

  for i from 1 to nops(objs) do
    if TrussMe:-IsRod(objs[i]) then
      # x coordinate of the joint location
      TrussMe:-Simplify(TrussMe:-ListPadding(eval(obj_coords[i]), 3))[1];
      # x coordinate of the joint location minus target length
      TrussMe:-Simplify(eval(%^2) - eval(objs[i]["length"]^~2));
      if  %% <> 0 and %% <> 0. and
          % <> 0 and % <> 0. then
        error "SUPPORT objects can only be applied at extremes of ROD objects"
      end if;
    end if;
    if TrussMe:-IsRod(objs[i]) and (constrained_dof[4..6] <> [0, 0, 0]) then
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
    if (S_stiffness(x) <> stiffness(x)) and m_WarningMode then
      WARNING("stiffness components not corresponding to constrained_dof are ignored");
    end if;
  else
    # Check for non zero stiffness on constrained dof
    if has(stiffness[remove(x -> x=0, ([seq(i, i = 1..6)]) *~ constrained_dof)], 0) then
      error "stiffness corresponding to constrained degrees of freedom cannot be zero";
    end if;
    # Check for zero stiffness on unconstrained dof
    S_stiffness := x -> stiffness *~ constrained_dof;
    if (S_stiffness(x) <> stiffness) and m_WarningMode then
      WARNING("stiffness components not corresponding to constrained_dof are ignored");
    end if;
  end if;

  S := table({
    "type"                = SUPPORT,
    "constrained_dof"     = constrained_dof,
    "admissible_loads"    = constrained_dof,
    "coordinates"         = [[0, 0, 0], op(map(TrussMe:-ListPadding, obj_coords, 3))],
    "name"                = name,
    "frame"               = RF,
    "targets"             = [m_earth["name"]] union TrussMe:-GetNames(objs),
    "variables"           = [],
    "forces"              = [],
    "moments"             = [],
    "constraint_loads"    = [],
    "support_reactions"   = [], # Expressed in support reference frame
    "stiffness"           = S_stiffness,
    "displacements"       = [],
    "frame_displacements" = LinearAlgebra:-IdentityMatrix(4)
    });

  # Build the temporary joint
  J_tmp := MakeJoint(name, constrained_dof, [m_earth, op(objs)], S["coordinates"], RF);

  S["variables"]        := J_tmp["variables"];
  S["forces"]           := J_tmp["forces"];
  S["moments"]          := J_tmp["moments"];
  S["constraint_loads"] := J_tmp["constraint_loads"];

  # Retrieve support force reactions
  sr_F_names := [FX, FY, FZ];
  for i from 1 to nops(S["forces"]) do
    if (S["forces"][i]["target"] = m_earth["name"]) then
      # Project forces in the support reference frame
      sr_F_values_tmp := TrussMe:-Project(S["forces"][i]["components"], ground, S["frame"]);
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
    if (S["moments"][i]["target"] = m_earth["name"]) then
      # Project moments in the support reference frame
      sr_M_values_tmp := TrussMe:-Project(S["moments"][i]["components"], ground, S["frame"]);
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

  return op(S);
end proc: # MakeSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsSupport := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Check if the object <obj> is a SUPPORT object.";

  return evalb(obj["type"] = SUPPORT);
end proc: # IsSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsCompliantSupport := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Check if the object <obj> is a SUPPORT object with compliant "
    "constraints.";

  local i;

  if TrussMe:-IsSupport(obj) then
    for i from 1 to 6 do
      if (obj["stiffness"](x)[i] <> infinity) and
         (obj["constrained_dof"][i] = 1) then
        return true;
      end if;
    end do;
  end if;
  return false;
end proc; # IsCompliantSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export CleanSupport := proc( # REVIEWED
  obj::SUPPORT,
  $)

  description "Clean SUPPORT object <obj> internal variables";

  obj["constraint_loads"]  := [];
  obj["support_reactions"] := [];
  obj["displacements"]     := [];
  return NULL;
end proc: # CleanSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeJoint := proc(
  name::string,                                               # Joint name
  constrained_dof::list,                                      # Constrained degree of freedom
  objs::list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH}), # Target objects
  coords::list,                                               # Joint locations
  RF::FRAME := ground,                                        # Reference frame
  {
    stiffness::{procedure, list(algebraic)} := [ # Stiffness components (default = infinite)
      infinity, infinity, infinity,
      infinity, infinity, infinity
    ] *~ constrained_dof,
    shell_objs:: # Objects to be considered connected to the shell of the joint
      list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH}) := [objs[1]]
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

  # Substitute -1 entries of coords with the corresponding object length
  obj_coords := [seq(`if`(coords[i] = -1, objs[i]["length"], coords[i]), i = 1..nops(coords))];

  for i from 1 to nops(objs) do

    if TrussMe:-IsRod(objs[i]) then
      # x coordinate of the joint location
      TrussMe:-ListPadding(eval(obj_coords[i]), 3)[1];
      # x coordinate of the joint location minus target length
      eval(eval(%^2) - eval(objs[i]["length"]^~2));
      if  %% <> 0 and %% <> 0. and
          % <> 0 and % <> 0. then
        error "JOINT objects can only be applied at extremes of ROD objects";
      end if;
    end if;
    if TrussMe:-IsRod(objs[i]) and (constrained_dof[4..6] <> [0, 0, 0]) then
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
    if (J_stiffness(x) <> stiffness(x)) and m_WarningMode then
      WARNING("stiffness components not corresponding to constrained_dof are ignored");
    end if;
  else
    # Check for non zero stiffness on constrained dof
    if has(stiffness[remove(x -> x = 0, ([seq(i, i = 1..6)]) *~ constrained_dof)], 0) then
      error "stiffness corresponding to constrained degrees of freedom cannot be zero";
    end if;
    # Check for zero stiffness on unconstrained dof
    J_stiffness := x -> stiffness *~ constrained_dof;
    if (J_stiffness(x) <> stiffness) and m_WarningMode then
      WARNING("stiffness components not corresponding to constrained_dof are ignored");
    end if;
  end if;

  J := table({
    "type"                = JOINT,
    "constrained_dof"     = constrained_dof,
    "admissible_loads"    = constrained_dof,
    "coordinates"         = map(TrussMe:-ListPadding, obj_coords, 3),
    "name"                = name,
    "frame"               = RF,
    "targets"             = TrussMe:-GetNames(objs),
    "shell_targets"       = TrussMe:-GetNames(shell_objs),
    "variables"           = [],
    "forces"              = [],
    "moments"             = [],
    "constraint_loads"    = [],
    "stiffness"           = J_stiffness,
    "displacements"       = [],
    "frame_displacements" = LinearAlgebra:-IdentityMatrix(4)
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
      TrussMe:-Project(jf_comp_cons, RF, objs[i]["frame"])
      .~ <op(objs[i]["admissible_loads"][1..3])>,
      list);
    # Check if there are reactions
    if (nops(jf_surv) <> 0) then
      # Create the reaction force between joint and obj
      # Force on obj
      JF_||(name)||_||(objs[i]["name"]) := TrussMe:-MakeForce(
        jf_comp_obj, obj_coords[i], objs[i], objs[i]["frame"]
        );
      # Force on joint
      JF_||(objs[i]["name"])||_||(name) := TrussMe:-MakeForce(
        -jf_comp_cons, 0, J, RF);
      # Use the non admissible loads to build the loads constraint
      constraint := convert(
        TrussMe:-Project(jf_comp_cons, RF, objs[i]["frame"])
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
      TrussMe:-Project(jm_comp_cons, RF, objs[i]["frame"])
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
        TrussMe:-Project(jm_comp_cons, RF, objs[i]["frame"])
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

  return op(J);
end proc: # MakeJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsJoint := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Check if the objecy <obj> is a JOINT object.";

  return evalb(obj["type"] = JOINT);
end proc: # IsJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsCompliantJoint := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Check if the object <obj> is a JOINT object with compliant "
    "constraints.";

  local i;

  if TrussMe:-IsJoint(obj) then
    for i from 1 to 6 do
      if (obj["stiffness"](x)[i] <> infinity) and
         (obj["constrained_dof"][i] = 1) then
        return true;
        break;
      end if;
    end do;
  end if;
  return false;
end proc; # IsCompliantJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export CleanJoint := proc( # REVIEWED
  obj::JOINT,
  $)

  description "Clean JOINT object <obj> internal variables.";

  obj["constraint_loads"] := [];
end proc: # CleanJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeRodPoints := proc( # REVIEWED
  rod_name::string,
  p_1::POINT,
  p_2::POINT,
  vec::Vector,
  {
    area::{algebraic, procedure} := infinity,
    material::MATERIAL           := TrussMe:-MakeMaterial()
  }, $)::ROD;

  description "Create a ROD object with inputs: object name <rod_name>, first "
    "point <p_1>, second point <p_2>, vector in XY-plane <vec>, optional "
    "section area <area> and material type <material>.";

  local ell, e_x, e_y, e_z;

  if (p_1 = p_2) then
    error "input points are the same";
  end if;

  # FIXME: does not work with symbolic vectors
  # if (TrussMe:-Norm2(vec) < 0) then
  #   error "input vector is null";
  # end if;

  ell := TrussMe:-Norm2(p_2 - p_1);
  e_x := (p_2 - p_1) /~ ell;
  e_y := convert(LinearAlgebra:-CrossProduct(<op(vec)>, <op(e_x)>), list);
  e_y := e_y /~ TrussMe:-Norm2(e_y);
  e_z := convert(LinearAlgebra:-CrossProduct(<op(e_x)>, <op(e_y)>), list);
  e_z := e_z /~ TrussMe:-Norm2(e_z);

  return TrussMe:-MakeRod(
    rod_name,
    ell,
    TrussMe:-Simplify(
      <<e_x[1], e_x[2], e_x[3], 0>|
       <e_y[1], e_y[2], e_y[3], 0>|
       <e_z[1], e_z[2], e_z[3], 0>|
       <p_1[1], p_1[2], p_1[3], 1>>
    ),
    parse("area")     = area,
    parse("material") = material
    );
end proc: # MakeRodPoints

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeRod := proc( # REVIEWED
  rod_name::string,
  ell::algebraic,
  RF::FRAME := ground,
  {
    area::{algebraic, procedure} := infinity,
    material::MATERIAL           := TrussMe:-MakeMaterial()
  }, $)::ROD;

  description "Create a ROD object with inputs: object name <rod_name>, reference "
    "length <ell>, optional reference frame <RF> in which the rod is defined, "
    "and optional section area <area> and material type <material>.";

  local area_proc;

  if type(area, procedure) then
    area_proc := area;
  else
    area_proc := x -> area;
  end if;

  return table({
    "type"                = ROD,
    "name"                = rod_name,
    "length"              = ell,
    "area"                = area_proc,
    "material"            = material,
    "frame"               = RF,
    "admissible_loads"    = [1, 0, 0, 0, 0, 0],
    "internal_actions"    = [],
    "displacements"       = [],
    "frame_displacements" = LinearAlgebra:-IdentityMatrix(4)
    });
end proc: # MakeRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsRod := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Check if the object <obj> is a ROD object.";

  return evalb(obj["type"] = ROD);
end proc: # IsRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export CleanRod := proc( # REVIEWED
  obj::ROD,
  $)

  description "Clean ROD object <obj> internal variables.";

  obj["internal_actions"] := [];
  obj["displacements"]    := [];
  return NULL;
end proc: # CleanRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeBeamPoints := proc( # REVIEWED
  beam_name::string,
  p_1::POINT,
  p_2::POINT,
  vec::Vector,
  {
    area::{algebraic, procedure}                   := infinity,
    timo_shear_coeff::{list(algebraic), procedure} := [5/6, 5/6],
    material::MATERIAL                             := TrussMe:-MakeMaterial(),
    I_xx::{algebraic, procedure}                   := infinity,
    I_yy::{algebraic, procedure}                   := infinity,
    I_zz::{algebraic, procedure}                   := infinity
  }, $)::BEAM;

  description "Create a BEAM object with inputs: object name <beam_name>, first "
    "point <p_1>, second point <p_2>, vector in XY-plane <vec>, optional "
    "section area <area>, optional Timoshenko shear coefficient <timo_shear_coeff>, "
    "optional material type <material>, optional section x-axis inertia <I_xx>, "
    "optional section y-axis inertia <I_yy> and optional section z-axis inertia "
    "<I_zz>.";

  local ell, e_x, e_y, e_z;

  if (p_1 = p_2) then
    error "input points are the same.";
  end if;

  # FIXME: does not work with symbolic vectors
  # if (TrussMe:-Norm2(vec) < 0) then
  #   error "input vector is null";
  # end if;

  ell := TrussMe:-Norm2(p_2 - p_1);
  e_x := (p_2 - p_1) /~ ell;
  e_y := convert(LinearAlgebra:-CrossProduct(<op(vec)>, <op(e_x)>), list);
  e_y := e_y /~ TrussMe:-Norm2(e_y);
  e_z := convert(LinearAlgebra:-CrossProduct(<op(e_x)>, <op(e_y)>), list);
  e_z := e_z /~ TrussMe:-Norm2(e_z);

  return TrussMe:-MakeBeam(
    beam_name,
    ell,
    TrussMe:-Simplify(
      <<e_x[1], e_x[2], e_x[3], 0>|
       <e_y[1], e_y[2], e_y[3], 0>|
       <e_z[1], e_z[2], e_z[3], 0>|
       <p_1[1], p_1[2], p_1[3], 1>>
    ),
    parse("area")             = area,
    parse("timo_shear_coeff") = timo_shear_coeff,
    parse("material")         = material,
    parse("I_xx")             = I_xx,
    parse("I_yy")             = I_yy,
    parse("I_zz")             = I_zz
  );
end proc: # MakeBeamPoints

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeBeam := proc( # REVIEWED
  beam_name::string,
  ell::algebraic,
  RF::FRAME := ground,
  {
    area::{algebraic, procedure}                   := infinity,
    timo_shear_coeff::{list(algebraic), procedure} := [5/6, 5/6],
    material::MATERIAL                             := TrussMe:-MakeMaterial(),
    I_xx::{algebraic, procedure}                   := infinity,
    I_yy::{algebraic, procedure}                   := infinity,
    I_zz::{algebraic, procedure}                   := infinity
  }, $)::BEAM;

  description "Create a BEAM object with inputs: object name <beam_name>, reference "
    "length <ell>, optional reference frame <RF> in which the rod is defined, and "
    "optional section area <area> and material type <material> and inertias on "
    "x- <I_xx>, y- <I_yy>, and z-axis <I_zz>.";

  local area_proc, timo_shear_coeff_proc, I_xx_proc, I_yy_proc, I_zz_proc;

  if type(area, procedure) then
    area_proc := area;
  else
    area_proc := x -> area;
  end if;

  if type(timo_shear_coeff, procedure) then
    timo_shear_coeff_proc := timo_shear_coeff;
  else
    timo_shear_coeff_proc := x -> timo_shear_coeff;
  end if;

  if type(I_xx, procedure) then
    I_xx_proc := I_xx;
  else
    I_xx_proc := x -> I_xx;
  end if;

  if type(I_yy, procedure) then
    I_yy_proc := I_yy;
  else
    I_yy_proc := x -> I_yy;
  end if;

  if type(I_zz, procedure) then
    I_zz_proc := I_zz;
  else
    I_zz_proc := x -> I_zz;
  end if;

  return table({
    "type"                = BEAM,
    "name"                = beam_name,
    "length"              = ell,
    "area"                = area_proc,
    "timo_shear_coeff"    = timo_shear_coeff_proc,
    "material"            = material,
    "inertias"            = [I_xx_proc, I_yy_proc, I_zz_proc],
    "frame"               = RF,
    "admissible_loads"    = [1, 1, 1, 1, 1, 1],
    "internal_actions"    = [],
    "displacements"       = [],
    "frame_displacements" = LinearAlgebra:-IdentityMatrix(4)
    });
end proc: # MakeBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsBeam := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Check if the object <obj> is a BEAM object.";

  return evalb(obj["type"] = BEAM);
end proc: # IsBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export CleanBeam := proc( # REVIEWED
  obj::BEAM,
  $)

  description "Clean BEAM object <obj> internal variables.";

  obj["internal_actions"] := [];
  obj["displacements"]    := [];
  return NULL;
end proc: # CleanBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeRigidBody := proc( # REVIEWED
  body_name::string,
  RF::FRAME := ground,
  {
    COM::list(algebraic) := [0, 0, 0],
    mass::algebraic      := 0
  },
  $)::RIGID_BODY;

  description "Create a RIGID_BODY object with inputs: object name <body_name>, "
    "reference frame <RF> in which the rigid body is defined, and optional "
    "center of mass position <COM> and mass <mass>.";

  return table({
    "type"                = RIGID_BODY,
    "name"                = body_name,
    "frame"               = RF,
    "COM"                 = COM,
    "mass"                = mass,
    "admissible_loads"    = [1, 1, 1, 1, 1, 1],
    "frame_displacements" = LinearAlgebra:-IdentityMatrix(4)
  });
end proc: # MakeRigidBody

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsRigidBody := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Check if the object <obj> is a RIGID_BODY object.";

  return evalb(obj["type"] = RIGID_BODY);
end proc: # IsRigidBody

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputeSpringDisplacement := proc(
  spring_load::algebraic,      # Load on the spring
  spring_stiffness::procedure, # Spring stiffness
  $)::algebraic;

  description "Compute the displacement of a spring give the load <spring_load> "
  "and spring stiffness <stiffness>.";

  local x, out, Dx;

  # Physics:-Assume(spring_load * Dx >= 0);

  # This works even for negative Dx
  out := RealDomain:-solve(
    spring_load = TrussMe:-Simplify(integrate(spring_stiffness(x), x = 0..Dx)), Dx, useassumptions = true
    );

  ##print("DISP EQ: ", spring_load = Simplify(integrate(spring_stiffness(x), x = 0..Dx)));

  #SolveTools:-Engine({
  #  spring_load = Simplify(integrate(spring_stiffness(x), x = 0..Dx))}, {Dx}, explicit):

  ##print("DISP SOL: ", %);
  #out := subs(Simplify(remove(x-> has(x,I), op~(%))), Dx); # FIXME deal with multiple solutions
  ##print("DISP SOL OUT: ", %);

  # Remove weird list notation inside a picewise function
  # FIXME: This works only if out is a single piece piecewise function
  if type(out, 'piecewise') then
    out := piecewise(op(map(x -> `if`(type(x, list), op(x), x), convert(out, list))));
  end if;

  return out;
end proc: # ComputeSpringDisplacement

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputeSpringEnergy := proc( # REVIEWED
  spring_load::algebraic,
  spring_stiffness::procedure,
  $)::algebraic;

  description "Compute the potential energy of a spring give the load <spring_load> "
  "and spring stiffness <stiffness>.";

  local disp;

  disp := TrussMe:-ComputeSpringDisplacement(spring_load, spring_stiffness);
  return TrussMe:-Simplify(
    integrate(integrate(spring_stiffness(x), x), x = 0..disp)
  );
end proc: # ComputeSpringEnergy

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputeSupportDisplacements := proc(
  obj::SUPPORT, # Support object
  $)

  description "Compute the displacements of the support <obj> from its "
    "support reactions.";

  local sup_disp, disp_vec, i, disp, x, sup_reac;

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

  if (m_VerboseMode > 0) then
    printf(
      "TrussMe:-ComputeSupportDisplacements(...): updating %s %s's displacements... ",
      obj["type"], obj["name"]
    );
  end if;

  obj["displacements"] := sup_disp;

  if (m_VerboseMode > 0) then
    printf("DONE\n");
  end if;

  TrussMe:-   out;
(procname);
end proc: # ComputeSupportDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputeJointDisplacements := proc( # REVIEWED
  obj::JOINT,
  sol::{list, set},
  $)

  description "Compute the displacements of the joint <obj> given the solution "
    "<sol>.";

  local jnt_disp, disp_vec, i, disp, x, jnt_load, f;

  jnt_disp := [];
  disp_vec := [tx, ty, tz, rx, ry, rz];
  jnt_load := [0, 0, 0, 0, 0, 0];

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

  if (m_VerboseMode > 0) then
    printf(
      "TrussMe:-ComputeJointDisplacements(...): updating %s %s's displacements... ",
      obj["type"], obj["name"]
    );
  end if;

  obj["displacements"] := jnt_disp;

  if (m_VerboseMode > 0) then
    printf("DONE\n");
  end if;

 return NULL;
end proc: # ComputeJointDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export MakeStructure := proc( # REVIEWED
  objs::{
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set( {BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  exts::{
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  } := [],
  {
    hyper_vars::{list ,set} := [],
    hyper_disp::{list, set} := [seq(0, 1..nops(hyper_vars))]
  }, $)::STRUCTURE;

  description "Create a STRUCTURE object with inputs: structure objects <objs>, "
    "external actions <exts>, optional hyperstatic variables <hyper_vars>, "
    "optional hyperstatic displacements <hyper_disp>.";

  local num_dof, names, candidate_hyp_vars, Graph, obj, S_ext;

  # Check for duplicate names
  names := [];
  for obj in objs do
    if member(obj["name"], names) then
      error "duplicate names found on structure objects";
    end if;
    names := names union [obj["name"]];
  end do;

  num_dof, Graph := TrussMe:-ComputeDOF(objs);

  if (num_dof < 0) then
    if (nops(hyper_vars) <> -num_dof) and m_WarningMode then
      candidate_hyp_vars := [];
      for obj in objs do
        if TrussMe:-IsSupport(obj) or TrussMe:-IsJoint(obj) then
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
      if (m_VerboseMode > 1) then
        printf(
          "TrussMe:-MakeStructure(...): hyperstatic structure detected with %d "
          "overconstrained directions.\n",
          abs(num_dof)
        );
      end if;
    end if;
  elif (num_dof > 0) and m_WarningMode then
    #error "not enough constraints in the structure";
    WARNING(
      "the structure is underconstrained with %1 unconstrained directions. "
      "Results computation may fail due to rigid body motions.",
      num_dof
    );
  else
    if (m_VerboseMode > 0) then
      printf("TrussMe:-MakeStructure(...): isostatic structure detected.\n");
    end if;
  end if;

  # Add gravity distributed load
  S_ext := exts;
  if (m_gravity <> [0, 0, 0]) then
    for obj in objs do
      if TrussMe:-IsRod(obj) and m_WarningMode then
        WARNING("TrussMe:-SolveStructure(...): gravity load is not supported for "
          "'ROD' type object %1", obj);
      elif TrussMe:-IsBeam(obj) then
        g_load||(obj["name"]) := TrussMe:-MakeQForce(
          (x -> m_gravity *~ obj["area"](x) *~ obj["material"]["density"]),
          obj,ground
          );
        S_ext := S_ext union {g_load||(obj["name"])};
      elif TrussMe:-IsRigidBody(obj) then
        g_load||(obj["name"]) := TrussMe:-MakeForce(
          m_gravity *~ obj["mass"], obj["COM"], obj, ground
          );
        S_ext := S_ext union {g_load||(obj["name"])};
      end if;
    end do;
  end if;

  return table({
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
end proc: # MakeStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsStructure := proc( # REVIEWED
  obj::anything,
  $)::boolean;

  description "Check if the object <obj> is a STRUCTURE object.";

  return type(obj, table) and evalb(obj["type"] = STRUCTURE); # FIXME adapt also for other types
end proc: # IsStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export CleanStructure := proc( # REVIEWED
  obj::STRUCTURE,
  $)

  description "Clean STRUCTURE object <obj> internal variables.";

  local i;

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
    if TrussMe:-IsBeam(obj[i]) then
      obj["objects"][i] := TrussMe:-CleanBeam(i);
    elif TrussMe:-IsRod(obj[i]) then
      obj["objects"][i] := TrussMe:-CleanRod(i);
    elif TrussMe:-IsSupport(obj[i]) then
      obj["objects"][i] := TrussMe:-CleanSupport(i);
    elif IsJoint(obj[i]) then
      obj["objects"][i] := TrussMe:-CleanJoint(i);
    end if;
  end do;
  return NULL;
end proc: # CleanStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputeDOF := proc( # REVIEWED
  objs::{
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set( {BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  }, $)::integer, function;

  description "Compute the degree of freedom of the input structure objects "
    "<objs>.";

  local dof, objs_tmp, obj, i, j, k, vertex, colors, G;

  dof      := 0;
  objs_tmp := objs union [m_earth];

  # Built connections graph
  vertex := [];
  colors := [];
  if (m_VerboseMode > 0) then
    printf("TrussMe:-ComputeDOF(...): checking structure connections... ");
  end if;
  for i from 1 to nops(objs_tmp) do
    vertex := vertex union [objs_tmp[i]["name"]];
    colors := colors union [ObjectColor(objs_tmp[i])];
  end do;
  G := GraphTheory:-Graph(vertex);
  GraphTheory:-HighlightVertex(G, vertex, colors);
  for i from 1 to nops(objs_tmp) do
    if TrussMe:-IsSupport(objs_tmp[i]) or TrussMe:-IsJoint(objs_tmp[i]) then
      for j from 1 to nops(objs_tmp) do
        if (member(objs_tmp[j]["name"], objs_tmp[i]["targets"])) then
          GraphTheory:-AddEdge(G, {objs_tmp[i]["name"], objs_tmp[j]["name"]});
        end if;
      end do;
    end if;
  end do;

  if (m_VerboseMode > 0) then
    printf("DONE\n");
    printf("TrussMe:-ComputeDOF(...): displaying connections graph... ");
    if (m_VerboseMode > 1) then
      print(GraphTheory:-DrawGraph(G), layout = tree);
    end if;
  end if;

  # Check graph connections
  if GraphTheory:-IsConnected(G) then
    if (m_VerboseMode > 0) then
      printf("DONE\n");
    end if;
  else
    print(GraphTheory:-DrawGraph(G), layout = tree);
    WARNING("unconnected elements detected in the structure");
  end if;

  if (m_VerboseMode > 0) then
    printf("TrussMe:-ComputeDOF(...): computing degrees of freedom... ");
  end if;

  for obj in objs_tmp do
    if TrussMe:-IsBeam(obj) then
      dof := dof + 6;
    elif TrussMe:-IsRigidBody(obj) then
      dof := dof + 6;
    elif TrussMe:-IsRod(obj) then
      dof := dof + 5;
    elif TrussMe:-IsJoint(obj) then
      dof := dof - add(
        obj["constrained_dof"][k], k = 1..6) * (nops(obj["targets"]) - 1
      );
    elif TrussMe:-IsSupport(obj) then
      dof := dof - add(
        obj["constrained_dof"][k], k = 1..6) * (nops(obj["targets"]) - 1
      );
    end if;
  end do;

  if (m_VerboseMode > 0) then
    printf("DONE (DOF = %d)\n", dof);
  end if;

  return dof, G;
end proc: # ComputeDOF

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DrawStructureGraph := proc( # REVIEWED
  obj::STRUCTURE,
  $)::function;

  description "Draw the connections graph of the STRUCTURE object <obj>.";

  return plots:-display(
    GraphTheory:-DrawGraph(obj["connections_graph"], layout = tree),
    title = "Structure connections graph"
  );
end proc: # DrawStructureGraph

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DrawStructureSparseMatrix := proc( # REVIEWED
  obj::STRUCTURE,
  {
    gauss_elimin::boolean := false
  },
  $)::procedure;

  description "Draw the sparse matrix for the equation system of STRUCTURE "
    "object <obj> and optionally apply Gaussian elimination to the matrix with "
    "the option <gauss_elimin>.";

  local A, B;

  A, B := LinearAlgebra:-GenerateMatrix(obj["equations"], obj["variables"]);
  if gauss_elimin then
    A := LinearAlgebra:-GaussianElimination(A);
  end if;

  return plots:-display(
    plot:-sparsematrixplot(A, matrixview),
    title = "Structure sparse matrix"
  );
end proc: # DrawStructureSparseMatrix

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export NewtonEuler := proc( # REVIEWED
  exts::{
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set({FORCE, MOMENT, QFORCE, QMOMENT})
  },
  obj::{BEAM, ROD, RIGID_BODY, SUPPORT, JOINT},
  {
    pole::POINT := [0, 0, 0],
    upper_lim::algebraic := obj["length"]
  }, $)::list;

  description "Compute the Newton-Euler static equilibrium equations given a set "
    "of external actions <exts>, and object to compute the equilibrium <obj>, "
    "the axial coordinate of the pole <pole>, and an optional upper limit of the "
    "integration <upper_lim>.";

  local eq_T, eq_R, i, x, arm, out;

  eq_T := [0, 0, 0];
  for i from 1 to nops(exts) do
    if exts[i]["target"] = obj["name"] then
      if TrussMe:-IsForce(exts[i]) then
        eq_T := eq_T + exts[i]["components"];
      elif TrussMe:-IsQForce(exts[i]) then
        eq_T := eq_T + convert(map(
          integrate, exts[i]["components"](x), x = 0..upper_lim
        ), list);
      end if;
    elif m_WarningMode then
      WARNING(
        "TrussMe:-NewtonEuler(...): %1 is not applied to %2.",
        exts[i], obj
      );
    end if;
  end do;

  eq_R := [0, 0, 0];
  for i from 1 to nops(exts) do
    if (exts[i]["target"] = obj["name"]) then
      if TrussMe:-IsMoment(exts[i]) then
        eq_R := eq_R + exts[i]["components"];
      elif TrussMe:-IsForce(exts[i]) then
        arm := <op(pole - exts[i]["coordinate"])>;
        eq_R := eq_R + convert(LinearAlgebra:-CrossProduct(
          <op(exts[i]["components"])>, arm
        ), list);
      elif TrussMe:-IsQForce(exts[i]) then
        arm := <op(pole - [x, 0, 0])>;
        eq_R := eq_R + map(integrate, convert(LinearAlgebra:-CrossProduct(
          <op(exts[i]["components"])(x)>, arm
        ), list), x = 0..upper_lim);
      elif TrussMe:-IsQMoment(FMQ[i]) then
        eq_R := eq_R + map(integrate, exts[i]["components"](x), x = 0..upper_lim);
      end if;
    elif m_WarningMode then
      WARNING(
        "TrussMe:-NewtonEuler(...): %1 is not applied to %2.",
        exts[i], obj
      );
    end if;
  end do;

  return eq_T union eq_R;
end proc: # NewtonEuler

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export SolveStructure := proc(
  struct::STRUCTURE, # Structure object
  {
    compute_internal_actions::boolean    := false, # Internal actions computation flag
    compute_displacements::boolean       := false, # Displacement computation flag
    compute_potential_energy::boolean    := false, # Potential energy computation flag
    compute_frame_displacements::boolean := false, # Frame displacements computation flag
    timoshenko_beam::boolean             := false, # Timoshenko beam flag
    implicit::boolean                    := false, # Implicit solution flag
    unveil_results::boolean              := true,  # Unveil results flag
    dummy_vars::{list, set}                := []     # Dummy variables
  }, $)::STRUCTURE;

  description "Solve the static equilibrium of a structure with inputs: "
    "structure <struct>, optional compute internal action enabling flag "
    "<compute_internal_actions>, optional compute displacement enabling flag "
    "<compute_displacements>, optional Timoshenko beam flag <timoshenko_beam>, "
    "optional implicit solution flag <implicit>, and optional unveil results "
    "flag <unveil_results>.";

  local g_load, S_obj, S_rigid, S_ext, S_support, S_joint, S_con_forces, vars,
    sol, obj, x, str_eq, str_vars, P_energy, veiling_idx, veils, i;

  # Clean structure
  TrussMe:-CleanStructure(struct);

  # Set veiling_label
  # FIXME veiling_idx := 1;
  # FIXME veiling_label := cat('_V', veiling_idx);

  # Parsing inputs
  S_obj        := {};
  S_rigid      := {};
  S_ext        := {};
  S_support    := {};
  S_joint      := {};
  S_con_forces := {};
  vars         := [];
  for obj in struct["objects"] do
    if TrussMe:-IsBeam(obj) or TrussMe:-IsRod(obj) then
      S_obj := S_obj union {eval(obj)};
    elif TrussMe:-IsRigidBody(obj) then
      S_rigid := S_rigid union {eval(obj)};
    elif TrussMe:-IsSupport(obj) then
      S_support    := S_support union {eval(obj)};
      S_con_forces := S_con_forces union obj["forces"] union obj["moments"];
      vars         := vars union obj["variables"];
    elif TrussMe:-IsJoint(obj) then
      S_joint      := S_joint union {eval(obj)};
      S_con_forces := S_con_forces union obj["forces"] union obj["moments"];
      vars         := vars union obj["variables"];
    end if;
    unassign('obj'); # FIXME: try to remove this line when the bug is fixed
  end do;

  S_ext := struct["external_actions"];

  # Set module local variable m_KeepVeiled
  m_KeepVeiled := not unveil_results;

  # Solve isostatic structure
  if (struct["dof"] >= 0) then

    if struct["dof"] > 0 and m_WarningMode then
      WARNING("TrussMe:-SolveStructure(...): structure is underconstrained. "
        "Trying to solve it anyway. Results computation may fail due to rigid "
        "body motions.");
    end if;

    if (m_VerboseMode > 0) then
      printf("TrussMe:-SolveStructure(...): solving the isostatic structure... ");
    end if;
    sol, str_eq, str_vars := TrussMe:-IsostaticSolver(
      S_obj union S_rigid union S_joint union S_support,
      S_ext union S_con_forces,
      vars,
      parse("implicit") = implicit
      );
    # Update Structure equations and variables
    struct["equations"] := str_eq;
    struct["variables"] := str_vars;
    if (m_VerboseMode > 0) then
      printf("DONE\n");
      if (m_VerboseMode > 1) then
        printf("TrussMe:-SolveStructure(...): updating support reactions fields... ");
      end if;
    end if;
    # Update support reactions properties
    for obj in S_support do
      obj["support_reactions"] := [
        seq(lhs(obj["support_reactions"][i]) = Subs(sol, rhs(obj["support_reactions"][i])),
        i = 1..nops(obj["support_reactions"]))
      ];
    end do;
    if (m_VerboseMode > 0) then
      printf("DONE\n");
    end if;

  # Solve hyperstatic structure
  elif (struct["dof"] < 0) then

    if (nops(struct["hyperstatic_variables"]) <> -struct["dof"]) then
      error "mismatch in the structure degrees of freedom, check the hyper"
        "static variables of the structure and update the structure object";
    end if;
    if (m_VerboseMode > 0) then
      printf("TrussMe:-SolveStructure(...): solving the hyperstatic structure... ");
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
    if (m_VerboseMode > 0) then
      printf("DONE\n");
      if (m_VerboseMode > 1) then
        printf("TrussMe:-SolveStructure(...): hyperstatic solver solution:\n");
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
    if (m_VerboseMode > 0) then
      printf("TrussMe:-SolveStructure(...): updating support reactions fields... ");
    end if;
    for obj in S_support do
    obj["support_reactions"] := Subs(sol, obj["support_reactions"]);
    end do;
    if (m_VerboseMode > 0) then
      printf("DONE\n");
    end if;
  end if;

  # Add veils
  if m_KeepVeiled and type(sol[-1], list) then
    struct["veils"] := struct["veils"] union sol[-1];
    # FIXME veiling_idx     := veiling_idx + 1;
    # FIXME veiling_label   := cat('_V', veiling_idx);
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
    if m_KeepVeiled then
      struct["veils"] := struct["veils"] union veils;
      # FIXME veiling_idx     := veiling_idx + 1;
      # FIXME veiling_label   := cat('_V', veiling_idx);
    end if;
  end if;

  return struct;
end proc: # SolveStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export HyperstaticSolver := proc( # REVIEWED
  objs::{
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  exts::{
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set({FORCE, MOMENT, QFORCE, QMOMENT})
  },
  vars::list,
  hyper_vars::list,
  hyper_disp::list,
  {
    timoshenko_beam::boolean := false,
    implicit::boolean        := false
  }, $)

  description "Solve hyperstatic structure with inputs objects <objs>, external "
    "actions <exts>, variables <vars>, hyperstatic variables <hyper_vars>, "
    "hyperstatic displacements <hyper_disp> and optional Timoshenko beam flag "
    "<timoshenko_beam>.";

  local hyper_eq, i, obj, iso_vars, iso_sol, iso_eq, hyper_sol, P_energy,
    S_objs, E_objs, sol;

  # Parse input objects and find objects with internal actions property
  S_objs := [seq(
    `if`(TrussMe:-IsBeam(objs[i]) or TrussMe:-IsRod(objs[i]), objs[i], NULL),
    i = 1..nops(objs))
    ];

  # Create a solution as function of the hyperstatic variables
  iso_vars := [seq(
    `if`(member(vars[i], hyper_vars), NULL, vars[i]),
    i = 1..nops(vars))
  ];
  iso_sol, iso_eq, iso_vars := TrussMe:-IsostaticSolver(
    objs, exts, iso_vars, parse("implicit") = implicit
  );

  # Compute internal actions
  ComputeInternalActions(S_objs, exts, iso_sol);

  # Extract the deformable objects
  E_objs := [seq(
    `if`(not TrussMe:-IsRigidBody(objs[i]), objs[i], NULL),
    i = 1..nops(objs))
  ];

  # Compute structure internal energy
  P_energy := TrussMe:-ComputePotentialEnergy(
    E_objs, iso_sol,
    parse("timoshenko_beam") = timoshenko_beam
  );
  P_energy := Simplify[60](P_energy);

  # Compute the hyperstatic equation
  if m_KeepVeiled and type(iso_sol[-1], list) then
    hyper_eq := TrussMe:-Diff~(
      P_energy, hyper_vars, parse("veils") = iso_sol[-1]
    ) =~ hyper_disp;
  else
    hyper_eq := diff~(P_energy, hyper_vars) =~ hyper_disp;
  end if;

  # Substitute Float(undefined) with 0 in the derivative of the potential energy
  # (this comes in case of non derivable piecewise functions in the compliant
  # joint stiffness)
  # FIXME: this is a temporary fix
  hyper_eq := TrussMe:-Simplify(eval(hyper_eq, Float(undefined) = 0));

  # Check for implicit solution flag
  if (implicit) then
    hyper_sol := iso_sol;
  else
    if (m_VerboseMode > 1) then
      printf("TrussMe:-HyperstaticSolver(...): solving the hyperstatic variables... ");
    end if;
    # Solve hyperstatic equations
    hyper_sol := op(RealDomain:-solve(hyper_eq, hyper_vars));
    # hyper_sol := op(SolveTools:-Engine({op(hyper_eq)}, {op(hyper_vars)}, explicit));
    # FIXME: if used must be adapted for non list equations (single equation)
    if (hyper_sol = NULL) then
      error "hyperstatic solution not found.";
    end if;

    # Substitute hyper_sol in P_energy
    P_energy := subs(hyper_sol, P_energy);

    if (m_VerboseMode > 1) then
      printf("DONE\n");
    end if;
    sol := hyper_sol union subs(hyper_sol, iso_sol);
  end if;

  if (_nresults = 4) then
    return sol, iso_eq union hyper_eq, iso_vars union hyper_vars, P_energy;
  else
    return sol;
  end
end proc: # HyperstaticSolver

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputePotentialEnergy := proc(
  objs::{
    list({BEAM, ROD, SUPPORT, JOINT}),
    set({BEAM, ROD, SUPPORT, JOINT})
  },
  sol::{list, set} := [],
  {
    timoshenko_beam::boolean := false,
    dummy_vars::{list, set}  := []
  }, $)

  description "Compute the internal potential energy of the structure given the "
    "objects <objs> and optional Timoshenko beam flag <timoshenko_beam>.";

  local dummy_vars_subs, obj, P, x, f, FJX, FJY, FJZ, MJX, MJY, MJZ, i;

  dummy_vars_subs := dummy_vars =~ [seq(0, i = 1..nops(dummy_vars))];

  P := 0;
  for obj in objs do
    if TrussMe:-IsBeam(obj) or TrussMe:-IsRod(obj) then
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
  return P;
end proc: # PotentialEnergy

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsostaticSolver := proc(
  objs::{
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  exts::{
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set({FORCE, MOMENT, QFORCE, QMOMENT})
  },
  vars::list,
  {
    implicit::boolean := false
  },
  $)

  description "Solve the isostatic structure equilibrium equation system given "
    "the structure objects <objs>, the external actions <exts> and the variables "
    "<vars> to solve.";

  local iso_eq, iso_eq_tmp, exts_comps, x, i, j, active_ext, iso_sol, A, B,
    iso_vars;

  # Compute structure equations
  if (m_VerboseMode > 1) then
    printf("TrussMe:-IsostaticSolver(...): computing the equilibrium equation "
      "for the isostatic structure...\n");
  end if;
  iso_eq := [];
  for i from 1 to nops(objs) do
    active_ext := {};
    for j from 1 to nops(exts) do
      if (exts[j]["target"] = objs[i]["name"]) then
        active_ext := active_ext union {exts[j]};
      end if;
    end do;
    iso_eq := iso_eq union TrussMe:-NewtonEuler(active_ext, objs[i]);
    # Add joints and supports constraint equations
    if TrussMe:-IsSupport(objs[i]) or TrussMe:-IsJoint(objs[i]) then
      iso_eq := iso_eq union objs[i]["constraint_loads"];
    end if;
  end do;

  # Remove NULL equations
  iso_eq := remove(x -> x = 0, TrussMe:-Simplify(iso_eq));

  # Remove equation related to rigid body motions
  # FIXME: for qloads should be different (this check can be skipped)
  iso_eq_tmp := iso_eq;
  exts_comps := map(x -> op(x["components"]), TrussMe:-GetObjsByType({FORCE, MOMENT}, exts));
  #iso_eq := remove(x -> (member(0., (abs~(exts_comps) -~ abs(x)) *~ 1.) and (not member(0., (abs~(vars) -~ abs(x)) *~ 1.))), iso_eq_tmp);
  iso_eq := remove(x -> (member(0., (abs~(exts_comps) -~ abs(x)) *~ 1.) and (not has(vars, indets(x)))), iso_eq_tmp);
  if m_WarningMode then
    if (nops(iso_eq_tmp) <> nops(iso_eq)) then
      WARNING("TrussMe:-IsostaticSolver(...): the following list of equations "
        "were removed because they are related to rigid body motions:\n%1",
        convert(iso_eq_tmp *~ 1., set) minus convert(iso_eq *~ 1., set));
    end if;
  end if;

  # Remove non used variables
  iso_vars := remove(x -> not has(iso_eq, x), vars);
  if m_WarningMode then
    if nops(iso_vars) <> nops(vars) then
      WARNING("TrussMe:-IsostaticSolver(...): the following list of variables "
        "were removed because they are not used in the structure equilibrium "
        "equations:\n%1", convert(vars, set) minus convert(iso_vars, set)
      );
    end if;
  end if;

  if (m_VerboseMode > 1) then
    printf("DONE\n");
    printf("TrussMe:-IsostaticSolver(...): structure equilibrium equations.\n");
    print(<op(iso_eq)>);
    printf("TrussMe:-IsostaticSolver(...): structure unknown variables.\n");
    print(iso_vars);
  end if;

  # Check for implicit solution flag
  if (implicit) then
    iso_sol := [];
  else
    # Matrix form
    A, B := LinearAlgebra:-GenerateMatrix(iso_eq, iso_vars);
    A := Matrix(A, storage = sparse);

    if nops(iso_eq) = nops(iso_vars) then

      if (m_VerboseMode > 1) then
        printf("TrussMe:-IsostaticSolver(...): matrix visualization of the linear system...\n");
        print(plots:-sparsematrixplot(A, matrixview));
      end if;
      if (m_VerboseMode > 1) then
        printf("TrussMe:-IsostaticSolver(...): computing the structure reaction forces... ");
      end if;
      # Solve structure equations (LinearSolver)
      #iso_sol := LinearSolver(iso_eq, iso_vars);
      iso_sol := op(RealDomain:-solve(iso_eq, iso_vars));
      # FIXME: this is a temporary fix (use LULEM when stable)
    else
      if m_WarningMode then
        WARNING("TrussMe:-IsostaticSolver(...): the system of equations is not "
          "consistent, trying to solve the system of equations anyway without "
          "LinearSolver");
      end if;
      # Solve structure equations (solve)
      iso_sol := op(RealDomain:-solve(iso_eq, iso_vars));
      # Append an empty list to the solution if the m_KeepVeiled flag
      # is set to true
      if m_KeepVeiled then
        iso_sol := iso_sol union [[]];
      end if;
    end if;

    if iso_sol = NULL then
      error "isostatic solution not found.";
    end if;

    if (m_VerboseMode > 1) then
      printf("DONE\n");
    end if;
  end if;

  if _nresults = 3 then
    return iso_sol, iso_eq, iso_vars;
  else
    return iso_sol;
  end
end proc: # IsostaticSolver

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputeInternalActions := proc(
  objs::{ # Structure objects
    list({BEAM, ROD}),
    set( {BEAM, ROD})
  },
  exts::{ # External actions
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  },
  sol::{list, set}, # Structure solution
  $)

  description "Programmatic computation of internal actions for structure"
    "objects with given external actions and structure solution.";

  local i, j, active_ext, subs_ext;

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
end proc: # ComputeInternalActions

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export InternalActions := proc(
  obj::{BEAM, ROD},
  exts::{
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set({FORCE, MOMENT, QFORCE, QMOMENT})
  },
  $)

  description "Programmatic computation of internal actions for structure "
    "object with given external actions, it returns the internal actions as "
    "function of the axial variable 'x'.";

  local i, ia, N_sol, Ty_sol, Tz_sol, Mx_sol, My_sol, Mz_sol, x, xx, xxx;


  # Clear assumptions if old assumptions
  # NOTE: assumptions help readability of the solution and improve computation
  # time, but results must be considered valid only in the assumed range
  Physics:-Assume(clear = x);

  # Meka some new assumptions
  Physics:-Assume(x > 0, x < obj["length"], x::real);

  # Compute internal actions for concentrated loads as effect overlay
  N_sol  := 0;
  Ty_sol := 0;
  Tz_sol := 0;
  Mx_sol := 0;
  My_sol := 0;
  Mz_sol := 0;
  for i from 1 to nops(exts) do
    if TrussMe:-IsForce(exts[i]) then
      N_sol  := N_sol  - TrussMe:-Simplify(piecewise(x >= exts[i]["coordinate"][1] and x <= obj["length"], exts[i]["components"][1]), piecewise);
      Ty_sol := Ty_sol + TrussMe:-Simplify(piecewise(x >= exts[i]["coordinate"][1] and x <= obj["length"], exts[i]["components"][2]), piecewise);
      Tz_sol := Tz_sol + TrussMe:-Simplify(piecewise(x >= exts[i]["coordinate"][1] and x <= obj["length"], exts[i]["components"][3]), piecewise);
      My_sol := My_sol - TrussMe:-Simplify(integrate(piecewise(x >= exts[i]["coordinate"][1] and x <= obj["length"], exts[i]["components"][3]), x = 0..x), piecewise);
      Mz_sol := Mz_sol + TrussMe:-Simplify(integrate(piecewise(x >= exts[i]["coordinate"][1] and x <= obj["length"], exts[i]["components"][2]), x = 0..x), piecewise);
    elif TrussMe:-IsMoment(exts[i]) then
      Mx_sol := Mx_sol - TrussMe:-Simplify(piecewise(x >= exts[i]["coordinate"][1] and x <= obj["length"], exts[i]["components"][1]), piecewise);
      My_sol := My_sol - TrussMe:-Simplify(piecewise(x >= exts[i]["coordinate"][1] and x <= obj["length"], exts[i]["components"][2]), piecewise);
      Mz_sol := Mz_sol - TrussMe:-Simplify(piecewise(x >= exts[i]["coordinate"][1] and x <= obj["length"], exts[i]["components"][3]), piecewise);
    elif TrussMe:-IsQForce(exts[i]) then
      N_sol  := N_sol  - TrussMe:-Simplify(integrate(exts[i]["components"](x)[1], x = 0..x), piecewise);
      Ty_sol := Ty_sol + TrussMe:-Simplify(integrate(exts[i]["components"](x)[2], x = 0..x), piecewise);
      Tz_sol := Tz_sol + TrussMe:-Simplify(integrate(exts[i]["components"](x)[3], x = 0..x), piecewise);
      My_sol := My_sol - TrussMe:-Simplify(integrate(integrate(exts[i]["components"](x)[3], x = 0..x), x = 0..x), piecewise);
      Mz_sol := Mz_sol + TrussMe:-Simplify(integrate(integrate(exts[i]["components"](x)[2], x = 0..x), x = 0..x), piecewise);
    elif TrussMe:-IsQMoment(FMQ[i]) then
      Mx_sol := Mx_sol - TrussMe:-Simplify(integrate(exts[i]["components"](x)[1], x = 0..x), piecewise);
      My_sol := My_sol - TrussMe:-Simplify(integrate(exts[i]["components"](x)[2], x = 0..x), piecewise);
      Mz_sol := Mz_sol - TrussMe:-Simplify(integrate(exts[i]["components"](x)[3], x = 0..x), piecewise);
    end if;
  end do;

  ia := [
    N  = unapply( N_sol, x),
    Ty = unapply(Ty_sol, x),
    Tz = unapply(Tz_sol, x),
    Mx = unapply(Mx_sol, x),
    My = unapply(My_sol, x),
    Mz = unapply(Mz_sol, x)
  ];

  if TrussMe:-IsRod(obj) then
    ia := [ia[1]];
  end if;

  if (m_VerboseMode > 1) then
  printf(
    "TrussMe:-InternalActions(...): updating %s %s's internal actions... ",
    obj["type"], obj["name"]
  );
  end if;

  obj["internal_actions"] := ia;

  if (m_VerboseMode > 1) then
    printf("DONE\n");
  end if;

  return NULL;
end proc: # InternalActions

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputeDisplacements := proc(
  objs::{ # Structural objects
    list({BEAM, ROD, SUPPORT, JOINT}),
    set({BEAM, ROD, SUPPORT, JOINT})
  },
  exts::{ # External loads
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set({FORCE, MOMENT, QFORCE, QMOMENT})
  },
  sol::{list, set}, # Solution of the structure
  {
    timoshenko_beam::boolean := false # Timoshenko beam flag
  }, $)

  description "Compute the structure displacements and rotations.";

  local obj, x, disp, rx_sol, ry_sol, rz_sol, ux_sol, uy_sol, uz_sol;

  # Cicle on the structure objects
  for obj in objs do
    # Beam
    if TrussMe:-IsBeam(obj) then
      Physics:-Assume(x > 0, x < obj["length"]);
      # Compute displacements
      rx_sol :=  integrate(subs(obj["internal_actions"](x), Mx(x)/(obj["material"]["shear_modulus"]*obj["inertias"][1](x))), x = 0..x);
      ry_sol :=  integrate(subs(obj["internal_actions"](x), My(x)/(obj["material"]["elastic_modulus"]*obj["inertias"][2](x))), x = 0..x);
      rz_sol :=  integrate(subs(obj["internal_actions"](x), Mz(x)/(obj["material"]["elastic_modulus"]*obj["inertias"][3](x))), x = 0..x);
      ux_sol :=  integrate(subs(obj["internal_actions"](x), N(x)/(obj["material"]["elastic_modulus"]*obj["area"](x))), x = 0..x);
      uy_sol :=  integrate(rz_sol, x = 0..x);
      uz_sol := -integrate(ry_sol, x = 0..x);
      if timoshenko_beam then
        uy_sol := uy_sol + integrate(subs(obj["internal_actions"](x), Ty(x)/(obj["timo_shear_coeff"](x)[1]*obj["material"]["shear_modulus"]*obj["area"](x))), x = 0..x);
        uz_sol := uz_sol + integrate(subs(obj["internal_actions"](x), Tz(x)/(obj["timo_shear_coeff"](x)[2]*obj["material"]["shear_modulus"]*obj["area"](x))), x = 0..x);
      end if;
      disp := [
        ux = unapply(ux_sol, x), uy = unapply(uy_sol, x), uz = unapply(uz_sol, x),
        rx = unapply(rx_sol, x), ry = unapply(ry_sol, x), rz = unapply(rz_sol, x)
      ];

      # Update object displacements
      obj["displacements"] := disp;

    # Rod
    elif TrussMe:-IsRod(obj) then
      Physics:-Assume(x > 0, x < obj["length"]);
      # Compute displacements
      ux_sol := integrate(subs(obj["internal_actions"](x), N(x)/(obj["material"]["elastic_modulus"]*obj["area"](x))), x = 0..x);
      disp := [
        ux = unapply(ux_sol, x)
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
end proc: # ComputeDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputePunctualDisplacement := proc(
  struct::STRUCTURE, # Structure
  objs::{ # Object on which the coordinates are defined
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  coords::list,     # Punctual coordinates defined in obj reference frame
  directions::list, # Displacement directions defined in obj reference frame
  RFs::list,        # Reference frames for the directions
  {
    timoshenko_beam::boolean := false, # Timoshenko beam flag
    unveil_results::boolean  := true   # Unveil results flag
  },
  $)::list, list;

  description "Compute the Structure <struct> punctual displacements of the "
    "object <obj> at the coordinates <coords> in the directions <directions>. "
    "The directions are defined in the reference frame <RFs>. Optional argument "
    "<timoshenko_beam> is a boolean flag to use Timoshenko beam model, <unveil_results> "
    "is a boolean flag to unveil the results.";

  local out, struct_copy, obj, objs_names, dummy_loads, subs_obj, obj_coords,
    obj_targets, x, subs_null_dummy, disp, i, j, d_coords, sw_tmp;

  # Substitute -1 entries of coords with the corresponding object length
  d_coords := [seq(`if`(coords[i] = -1, objs[i]["length"], coords[i]), i = 1..nops(coords))];

  # Set module local variable m_KeepVeiled
  m_KeepVeiled := not unveil_results;

  # Create a copy of the structure
  struct_copy := TrussMe:-CopyStructure(struct);

  # Get objects names
  objs_names := TrussMe:-GetNames(objs);

  # Disable warnings temporarily FIXME: this is a temporary fix
  sw_tmp := m_WarningMode;
  m_WarningMode := false;

  # Replace Rods with Beams to be able to compute the displacements in all directions
  for obj in map(TrussMe:-GetObjByName, objs_names, struct_copy["objects"]) do
    if TrussMe:-IsRod(obj) then
      subs_obj := TrussMe:-MakeBeam(obj["name"], obj["length"], obj["frame"], parse("area") = obj["area"], parse("material") = obj["material"]);
      # Remove load on unconstrained direction
      subs_obj["admissible_loads"] := [1, 1, 1, 0, 1, 1];
      # Replace object in struct_copy
      struct_copy["objects"] := remove(x -> x["name"] = obj["name"], struct_copy["objects"]);
      struct_copy["objects"] := struct_copy["objects"] union {eval(subs_obj)};
    end if;
  end do;

  # Update struct_copy supports and joints for the new objects
  for obj in struct_copy["objects"] do
    if TrussMe:-IsSupport(obj) then
      # Re-Make support to generate new loads and constraint compliant with substituted objects (joint is made because earth is already in the list of targets)
      obj_targets := map(TrussMe:-GetObjByName, remove(x -> x = m_earth["name"], obj["targets"]), struct_copy["objects"]);
      obj_coords  := obj["coordinates"][2..-1];
      subs_obj := MakeSupport(obj["name"], obj["constrained_dof"], obj_targets, obj_coords, obj["frame"], parse("stiffness") = obj["stiffness"]);
      # Replace object in struct_copy
      struct_copy["objects"] := remove(x -> x["name"] = obj["name"], struct_copy["objects"]);
      struct_copy["objects"] := struct_copy["objects"] union {eval(subs_obj)};
    elif TrussMe:-IsJoint(obj) then
      # Re-Make joint to generate new loads and constraint compliant with substituted objects
      obj_targets := map(TrussMe:-GetObjByName, obj["targets"], struct_copy["objects"]);
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
    `if`(directions[i,1] = 1, TrussMe:-MakeForce( [dFx_||i,0,0], d_coords[i], TrussMe:-GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL),
    `if`(directions[i,2] = 1, TrussMe:-MakeForce( [0,dFy_||i,0], d_coords[i], TrussMe:-GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL),
    `if`(directions[i,3] = 1, TrussMe:-MakeForce( [0,0,dFz_||i], d_coords[i], TrussMe:-GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL),
    `if`(directions[i,4] = 1, TrussMe:-MakeMoment([dMx_||i,0,0], d_coords[i], TrussMe:-GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL),
    `if`(directions[i,5] = 1, TrussMe:-MakeMoment([0,dMy_||i,0], d_coords[i], TrussMe:-GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL),
    `if`(directions[i,6] = 1, TrussMe:-MakeMoment([0,0,dMz_||i], d_coords[i], TrussMe:-GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL)
    ]);

    # Null dummy loads substitution list
    subs_null_dummy := subs_null_dummy union ([dFx_||i, dFy_||i, dFz_||i, dMx_||i, dMy_||i, dMz_||i] =~ [0, 0, 0, 0, 0, 0]);

    # Add dummy loads to the structure copy
    struct_copy["external_actions"] := struct_copy["external_actions"] union dummy_loads;
  end do;

  m_StoredData := m_StoredData union subs_null_dummy;

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
  m_WarningMode := sw_tmp;

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

   struct_copy["veils"] := subs(subs_null_dummy, struct_copy["veils"]);

  # Simplify output
  out := TrussMe:-Simplify(out);

  if _nresults = 1 then
    return out;
  else
    return out, struct_copy["veils"];
  end if;
end proc: # ComputePunctualDisplacement

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputeObjectFrameDisplacements := proc(
  struct::STRUCTURE, # Structure to compute the total displacements
  {
    timoshenko_beam::boolean := false, # Timoshenko beam flag
    unveil_results::boolean  := true   # Unveil results flag
  },
  $)::list;

description "Compute the total displacements of the structure <struct>.";

  local i, RF_nt, nx, ny, nz, theta, subs_n, subs_t, disp, veils, objs;

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

  return veils;
end proc; # ComputeObjectFrameDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export CopyStructure := proc(
  struct::STRUCTURE, # Structure to copy
  $)::STRUCTURE;

description "Create a copy of the structure <struct> and its objects.";

  local struct_copy, obj, action;

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

  return struct_copy;
end proc: # CopyStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export LinearSolver := proc(
  eqns::{list, set}, # Equations
  vars::{list, set}, # Variables
  $)::list;

  description "Solve the linear system of equations <eqns> for the variables "
  "<vars>.";

  local T, sol, sol_tmp, A, B, _Q, PivotingStrategy;

  # Matrix form of the linear system
  A, B := LinearAlgebra:-GenerateMatrix(eqns, vars);

  if (m_VerboseMode > 0) then
    printf("TrussMe:-LinearSolver(...): performing LU decomposition... ");
  end if;

  # LU decomposition
  T := LU(A);

  if (m_VerboseMode > 0) then
    printf("DONE (rank = %d)\n", T["rank"]);
    printf("TrussMe:-LinearSolver(...): solved linear system... ");
  end if;

  # Solve linear system
  sol_tmp := SolveLinearSystem(T, b);

  if (m_VerboseMode > 0) then
    printf("DONE\n");
    printf("TrussMe:-LinearSolver(...): substituting veils... ");
  end if;

  # Substitute veils to solution
  if m_KeepVeiled then
    # Remove indexed type from veils
    m_LEM:-VeilList();
    lhs~(%) =~ map2(op, 0, lhs~(%)) ||~ __ ||~ (op~(lhs~(%)));
    # Substitutution
    sol := convert(vars =~ subs(%, sol_tmp), list) union [subs(%, %%)];
  else
    sol := convert(vars =~ m_LEM:-VeilSubs(sol_tmp), list);
  end if;
  m_LEM:-VeilForget();

  if (m_VerboseMode > 0) then
    printf("DONE\n");
  end if;

  return TrussMe:-Simplify(sol);
end proc: # LinearSolver

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ObjectColor := proc( # REVIEWED
  obj::{BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH},
  $)::string;

  description "Return the color of the object <obj>.";

  if TrussMe:-IsBeam(obj) then
    return m_BeamColor;
  elif TrussMe:-IsRod(obj) then
    return m_RodColor;
  elif TrussMe:-IsRigidBody(obj) then
    return m_RigidBodyColor;
  elif TrussMe:-IsCompliantSupport(obj) then
    return m_CompliantSupportColor;
  elif TrussMe:-IsSupport(obj) then
    return m_SupportColor;
  elif TrussMe:-IsCompliantJoint(obj) then
    return m_CompliantJointColor;
  elif TrussMe:-IsJoint(obj) then
    return m_JointColor;
  elif TrussMe:-IsEarth(obj) then
    return m_EarthColor;
  else
    error "invalid object detected.";
  end if;
end proc: # ObjectColor

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotRigidBody := proc( # REVIEWED
  obj::RIGID_BODY,
  joints::{
    list({SUPPORT, JOINT}),
    set({SUPPORT, JOINT})
  },
  c_loads::{
    list({FORCE, MOMENT}),
    set({FORCE, MOMENT})
  },
  {
    data::{list(`=`), set(`=`)} := []
  },
  $)::function;

  description "Plot the RIGID_BODY object <obj>.";

  local p_1, p_2, i, idx, lines, load;

  lines := [];
  p_1 := subs(op(data),
    TrussMe:-Project([op(obj["COM"]), 1], obj["frame"], ground)
  );
  for i in joints do
    member(obj["name"], i["targets"], 'idx');
    p_2 := subs(op(data),
      TrussMe:-Project([op(i["coordinates"][idx]), 1], obj["frame"], ground)
    );
    lines := lines union [
      plottools:-line(
        convert(p_1[1..3], list),
        convert(p_2[1..3], list),
        thickness = 6
    )];
  end do;

  for load in c_loads do
    p_2 :=  subs(op(data),
      TrussMe:-Project([op(load["coordinate"]), 1], obj["frame"], ground)
    );
    lines := lines union [
      plottools:-line(
        convert(p_1[1..3], list),
        convert(p_2[1..3], list),
        thickness = 6
    )];
  end do;

  return plots:-display(
    lines,
    linestyle = solid,
    color     = TrussMe:-ObjectColor(obj),
    scaling   = constrained
  );
end proc: # PlotRigidBody

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotDeformedRigidBody := proc(
  obj::RIGID_BODY, # Rigid body to be plot
  joints::{        # Joint and support objects
    list({SUPPORT, JOINT}),
    set({SUPPORT, JOINT})
  },
  c_loads::{ # Concentrated loads
    list({FORCE, MOMENT}),
    set({FORCE, MOMENT})
  },
  {
    data::{list(`=`), set(`=`)} := [], # Substitutions
    scaling::numeric            := 1.0 # Scaling factor
  },
  $)::function;

  description "Plot the deformed RIGID_BODY object <obj>.";

  local P1, P2, rfd, js, idx, lines, load, out;

  lines := [];
  rfd := subs(op(data), obj["frame"] . ((obj["frame_displacements"] - LinearAlgebra:-IdentityMatrix(4)) *~ scaling + LinearAlgebra:-IdentityMatrix(4)));
  P1 := subs(op(data), TrussMe:-Project([op(obj["COM"]), 1], rfd, ground));
  for js in joints do
    member(obj["name"], js["targets"], 'idx');
    P2 := subs(op(data), TrussMe:-Project([op(js["coordinates"][idx]), 1], rfd, ground));
    lines := lines union [plottools:-line(convert(P1[1..3], list), convert(P2[1..3], list), thickness = 6)];
  end do;

  for load in c_loads do
    P2 :=  subs(op(data), TrussMe:-Project([op(load["coordinate"]), 1], rfd, ground));
    lines := lines union [plottools:-line(convert(P1[1..3], list), convert(P2[1..3], list), thickness = 6)];
  end do;

  out := plots:-display(lines, linestyle = solid, color = ObjectColor(obj), parse("scaling") = constrained);

  return out;
end proc: # PlotDeformedRigidBody

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotBeam := proc( # REVIEWED
  obj::BEAM,
  {
    data::{list(`=`), set(`=`)} := []
  },
  $)::function;

  description "Plot the BEAM object <obj>.";

  local p_1, p_2;

  p_1 := subs(op(data), TrussMe:-Origin(obj["frame"]));
  p_2 := subs(op(data),
    TrussMe:-Project([obj["length"], 0, 0, 1], obj["frame"], ground)
  );

  return plots:-display(
    plottools:-line(
      convert(p_1[1..3], list),
      convert(p_2[1..3], list),
      thickness = 6
    ),
    linestyle = solid,
    color     = TrussMe:-ObjectColor(obj),
    scaling   = constrained
  );
end proc: # PlotBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotDeformedBeam := proc( # REVIEWED
  obj::BEAM,
  {
    data::{list(`=`), set(`=`)} := [],
    scaling::numeric            := 1.0
  },
  $)::function;

  description "Plot the BEAM object <obj>.";

  local sc, rfd, x;

  rfd := subs(op(data), obj["frame"] . ((obj["frame_displacements"] -
    LinearAlgebra:-IdentityMatrix(4)) *~ scaling + LinearAlgebra:-IdentityMatrix(4)));

  sc := subs(op(data), TrussMe:-Project(
    subs(obj["displacements"](x),
      [x, 0, 0, 0] +~ [ux(x) *~ scaling, uy(x) *~ scaling, uz(x) *~ scaling, 1]
      ), rfd, ground)[1..3]);

  return plots:-display(
    plots:-spacecurve(
      sc, x = subs(op(data), 0..obj["length"]),
      thickness = 6
    ),
    linestyle        = solid,
    color            = TrussMe:-ObjectColor(obj),
    parse("scaling") = constrained
  );
end proc: # PlotDeformedBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotRod := proc( # REVIEWED
  obj::ROD,
  {
    data::{list(`=`), set(`=`)} := []
  },
  $)::function;

  description "Plot the ROD object <obj> given a list or set of substitutions "
    "data <data>.";

  local p_1, p_2, out;

  p_1 := subs(op(data), TrussMe:-Origin(obj["frame"]));
  p_2 := subs(op(data),
    TrussMe:-Project([obj["length"], 0, 0, 1], obj["frame"], ground)
  );

  return plots:-display(
    plottools:-line(
      convert(p_1[1..3], list),
      convert(p_2[1..3], list),
      thickness = 4
    ),
    linestyle = solid,
    color     = TrussMe:-ObjectColor(obj),
    scaling   = constrained
  );
end proc: # PlotRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotDeformedRod := proc(
  obj::ROD,
  {
    data::{list(`=`), set(`=`)} := [],
    scaling::numeric            := 1
  },
  $)::function;

  description "Plot the ROD object <obj> give the list of <targets> and "
    "a list or set of substitutions data <data>.";

  local P1, P2, rfd;

  rfd := obj["frame"] . ((obj["frame_displacements"] -
    LinearAlgebra:-IdentityMatrix(4)) *~ scaling + LinearAlgebra:-IdentityMatrix(4));

  P1 := subs(op(data), TrussMe:-Origin(rfd));
  P2 := subs(op(data), TrussMe:-Project(
    [obj["length"] + subs(obj["displacements"](obj["length"]
    ), ux(obj["length"]) *~ scaling), 0, 0, 1], rfd, ground));

  return plots:-display(
    plottools:-line(
      convert(P1[1..3], list),
      convert(P2[1..3], list),
      thickness = 4
    ),
    linestyle        = solid,
    color            = TrussMe:-ObjectColor(obj),
    parse("scaling") = constrained
  );
end proc: # PlotDeformedRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotJoint := proc( # REVIEWED
  obj::JOINT,
  targets::{
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  {
    data::{list(`=`), set(`=`)} := []
  },
  $)::function;

  description "Plot the JOINT object <obj> give the list of <targets> and "
    "a list or set of substitutions data <data>.";

  local O;

  O := subs(op(data), TrussMe:-Origin(
    TrussMe:-GetObjByName(obj["targets"][1], targets)["frame"].
    TrussMe:-Translate(op(TrussMe:-ListPadding(obj["coordinates"][1],3)))
    ));

  return plots:-display(
    plottools:-point(
      convert(O[1..3], list),
      symbol     = solidsphere,
      symbolsize = 20
    ),
    linestyle = solid,
    color     = TrussMe:-ObjectColor(obj),
    scaling   = constrained
  );
end proc: # PlotJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotDeformedJoint := proc( # REVIEWED
  obj::JOINT,
  targets::{
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  {
    data::{list(`=`), set(`=`)} := [],
    scaling::numeric            := 1
  },
  $)::function;

  description "Plot the JOINT object <obj>.";

  local O, rfd;

  #TODO: add compliant joint deformation

  rfd := ((obj["frame_displacements"] - LinearAlgebra:-IdentityMatrix(4)) *~ scaling + LinearAlgebra:-IdentityMatrix(4));

  # FIXME: joint target may be a joint itself
  O := subs(op(data), TrussMe:-Origin(
     TrussMe:-GetObjByName(obj["targets"][1], targets)["frame"].
     TrussMe:-Translate(op(TrussMe:-ListPadding(obj["coordinates"][1],3))))[1..3] +~
     TrussMe:-Project(TrussMe:-Origin(rfd)[1..3], obj["frame"], ground)
    );

  return plots:-display(
    plottools:-point(
      convert(O, list),
      symbol     = solidsphere,
      symbolsize = 20
    ),
    linestyle        = solid,
    color            = TrussMe:-ObjectColor(obj),
    parse("scaling") = constrained);
end proc: # PlotDeformedJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotSupport := proc( # REVIEWED
  obj::SUPPORT,
  targets::{
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  {
    data::{list(`=`), set(`=`)} := []
  },
  $)::function;

  description "Plot the SUPPORT object <obj>.";

  local O, out;

  if (nops(obj["targets"]) > 1) then
    O := subs(op(data), TrussMe:-Origin(
      TrussMe:-GetObjByName(obj["targets"][2], targets)["frame"].
      TrussMe:-Translate(op(TrussMe:-ListPadding(obj["coordinates"][2],3)))
      ));
  else
    O := subs(op(data), TrussMe:-Origin(
      m_earth["frame"].
      TrussMe:-Translate(op(TrussMe:-ListPadding(obj["coordinates"][1],3)))
      ));
  end if;

  return plots:-display(
    plottools:-point(
      convert(O[1..3], list),
      symbol     = solidbox,
      symbolsize = 20
    ),
    linestyle = solid,
    color     = TrussMe:-ObjectColor(obj),
    scaling   = constrained
  );
end proc: # PlotSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotDeformedSupport := proc( # REVIEWED
  obj::SUPPORT,
  targets::{
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  {
    data::{list(`=`), set(`=`)} := [],
    scaling::numeric            := 1.0
  },
  $)::function;

  description "Plot the deformed SUPPORT object <obj>.";

  local O, out;

  # TODO: add compliant support deformation

  if (nops(obj["targets"]) > 1) then
    O := subs(op(data), TrussMe:-Origin(
      TrussMe:-GetObjByName(obj["targets"][2], targets)["frame"].
      TrussMe:-Translate(op(TrussMe:-ListPadding(obj["coordinates"][2],3)))
      ));
  else
    O := subs(op(data), TrussMe:-Origin(
      m_earth["frame"].
      TrussMe:-Translate(op(TrussMe:-ListPadding(obj["coordinates"][1],3)))
      ));
  end if;

  return plots:-display(
    plottools:-point(
      convert(O[1..3], list),
      symbol     = solidbox,
      symbolsize = 20
    ),
    linestyle        = solid,
    color            = TrussMe:-ObjectColor(obj),
    parse("scaling") = constrained
  );
end proc: # PlotDeformedSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotStructure := proc( # REVIEWED
  str::STRUCTURE,
  {
    data::{list(`=`), set(`=`)} := []
  },
  $)::{function, list(function)};

  description "Plot the STRUCTURE object <str> given a list of substitutions "
    "<data>.";

  local disp, rb_joints, rb_loads, obj;

  disp := []:
  for obj in str["objects"] do
    if TrussMe:-IsBeam(obj) then
      disp := disp union [
        TrussMe:-PlotBeam(obj, parse("data") = data)
      ];
    elif TrussMe:-IsRod(obj) then
      disp := disp union [
        TrussMe:-PlotRod(obj, parse("data") = data)
      ];
    elif TrussMe:-IsSupport(obj) then
      map(GetObjByName, obj["targets"], str["objects"]);
      disp := disp union [TrussMe:-PlotSupport(obj, %, parse("data") = data)];
    elif TrussMe:-IsJoint(obj) then
      map(GetObjByName, obj["targets"], str["objects"]);
      disp := disp union [TrussMe:-PlotJoint(obj, %, parse("data") = data)];
    elif TrussMe:-IsRigidBody(obj) then
      TrussMe:-GetObjsByType(['JOINT', 'SUPPORT'], str["objects"]);
      rb_joints := remove(x -> (not member(obj["name"], x["targets"])), %);
      TrussMe:-GetObjsByType(['FORCE', 'MOMENT'], str["external_actions"]);
      rb_loads := remove(x -> obj["name"] <> x["target"], %);
      disp := disp union [
        TrussMe:-PlotRigidBody(obj, rb_joints, rb_loads, parse("data") = data)
      ];
    end if;
  end do;

  return plots:-display(
    disp,
    axes    = boxed,
    scaling = constrained,
    labels  = ['x', 'y', 'z']
  );
end proc: # PlotStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotDeformedStructure := proc(
  str::STRUCTURE,
  {
    data::{list(`=`), set(`=`)} := [],
    scaling::numeric            := 1
  },
  $)::{function, list(function)};

  description "Plot the deformed STRUCTURE object <str> given a list of "
    "substitutions <data> and a scaling factor <scaling>.";

  local disp, rb_joints, rb_loads, obj;

  # Check if displacements and frame displacements are solved
  if not str["displacements_solved"] or not str["frame_displacements_solved"] then
    error "displacements and frame displacements must be solved before plotting "
      "the deformed structure.";
  end if;

  disp := []:
  for obj in str["objects"] do
    if TrussMe:-IsBeam(obj) then
      disp := disp union [
        TrussMe:-PlotDeformedBeam(obj, parse("data") = data, parse("scaling") = scaling)
      ];
    elif TrussMe:-IsRod(obj) then
      disp := disp union [
        TrussMe:-PlotDeformedRod(obj, parse("data") = data, parse("scaling") = scaling)
      ];
    elif TrussMe:-IsSupport(obj) then
      disp := disp union [
        TrussMe:-PlotDeformedSupport(obj, map(TrussMe:-GetObjByName, obj["targets"], str["objects"]), parse("data") = data, parse("scaling") = scaling)
      ];
    elif TrussMe:-IsJoint(obj) then
      disp := disp union [
        TrussMe:-PlotDeformedJoint(obj, map(TrussMe:-GetObjByName, obj["targets"], str["objects"]), parse("data") = data, parse("scaling") = scaling)
      ];
    elif TrussMe:-IsRigidBody(obj) then
      TrussMe:-GetObjsByType(['JOINT', 'SUPPORT'], str["objects"]);
      rb_joints := remove(x -> (not member(obj["name"], x["targets"])), %);
      TrussMe:-GetObjsByType(['FORCE', 'MOMENT'], str["external_actions"]);
      rb_loads := remove(x -> obj["name"] <> x["target"], %);
      disp := disp union [
        TrussMe:-PlotDeformedRigidBody(obj, rb_joints, rb_loads, parse("data") = data, parse("scaling") = scaling)
      ];
    end if;
  end do;

  return plots:-display(
    disp, axes = boxed, parse("scaling") = constrained, labels = ['x', 'y', 'z']
  );
end proc: # PlotDeformedStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsInsideJoint := proc( # REVIEWED
  obj::JOINT,
  pnt::POINT,
  tol::numeric := 1e-4,
  $)::boolean;

  description "Check if the point <pnt> is inside the JOINT <obj> within a "
    "tolerance <tol>.";

  local O;

  if not (nops(pnt) = 3) then
    error "the input point must be a list of 3 elements.";
  end if;

  if (nops(obj["targets"]) > 1) then
    O := TrussMe:-Origin(
      TrussMe:-GetObjByName(obj["targets"][1], targets)["frame"]. # FIX 'targets' as a symbol???
      TrussMe:-Translate(obj["coordinates"][1], 0, 0)
      );
  elif m_WarningMode then
    WARNING("The support has no targets");
  end if;

  return evalb(norm(pnt - O) <= tol);
end proc: # IsInsideJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsInsideSupport := proc( # REVIEWED
  obj::SUPPORT,
  pnt::POINT,
  tol::numeric := 1e-4,
  $)::boolean;

  description "Check if the point <p> is inside the SUPPORT <obj> within a "
    "tolerance <tol>.";

  local O;

  if not (nops(pnt) = 3) then
    error "the input point must be a list of 3 elements.";
  end if;

  if (nops(obj["targets"]) > 1) then
    O := TrussMe:-Origin(
      TrussMe:-GetObjByName(obj["targets"][2], targets)["frame"]. # FIX 'targets' as a symbol???
      TrussMe:-Translate(obj["coordinates"][2], 0, 0)
      );
  elif m_WarningMode then
    WARNING("the support has no targets");
  end if;

  return evalb(norm(pnt - O) <= tol);
end proc: # IsInsideSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsInsideRod := proc( # REVIEWED
  obj::ROD,
  pnt::POINT,
  $)::boolean;

  description "Check if the point <pnt> is inside the ROD <obj>.";

  local O, V, W;

  if not (nops(pnt) = 3) then
    error "the input point must be a list of 3 elements.";
  end if;

  O := TrussMe:-Origin(obj["frame"]);
  V := obj["frame"].TrussMe:-Translate(obj["length"], 0, 0) - O;
  W := pnt - O;
  return evalb(dot(W, V) >= 0) and (dot(W, V) <= dot(V, V));
end proc: # IsInsideRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsInsideBeam := proc( # REVIEWED
  obj::BEAM,
  pnt::POINT,
  $)::boolean;

  description "Check if the point <pnt> is inside the BEAM <obj>.";

  local O, V, W, out;

  if not (nops(pnt) = 3) then
    error "the input point must be a list of 3 element.";
  end if;

  O := TrussMe:-Origin(obj["frame"]);
  V := obj["frame"].TrussMe:-Translate(obj["length"], 0, 0) - O;
  W := pnt - O;
  return evalb((dot(W, V) >= 0) and (dot(W, V) <= dot(V, V)));
end proc: # IsInsideBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsInsideStructure := proc( # REVIEWED
  obj::STRUCTURE,
  pnt::POINT,
  $)::boolean;

  description "Check if the point <pnt> is inside the STRUCTURE <obj>.";

  local i, out;

  out := false;
  for i in str["objects"] do
    if TrussMe:-IsJoint(i) then
      out := out or TrussMe:-IsInsideJoint(i, pnt);
    elif TrussMe:-IsSupport(i) then
      out := out or TrussMe:-IsInsideSupport(i, pnt);
    elif TrussMe:-IsBeam(i) then
      out := out or TrussMe:-IsInsideBeam(i, pnt);
    elif TrussMe:-IsRod(i) then
      out := out or TrussMe:-IsInsideRod(i, pnt);
    else
      error "unknown object type.";
    end if;
    if out then
      break;
    end if;
  end do;
  return out;
end proc: # IsInsideStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module:

# That's all folks!
