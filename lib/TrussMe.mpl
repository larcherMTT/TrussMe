# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                     _____                   __  __                          #
#                    |_   _| __ _   _ ___ ___|  \/  | ___                     #
#                      | || '__| | | / __/ __| |\/| |/ _ \                    #
#                      | || |  | |_| \__ \__ \ |  | |  __/                    #
#                      |_||_|   \__,_|___/___/_|  |_|\___|                    #
#                A Maple Library for Truss Elements Structures                #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Current version authors:
#   Matteo Larcher (University of Trento)
#   Davide Stocco (University of Trento)
#
# License: BSD 3-Clause License
#
# This is a module for the 'TrussMe' (A Maple Library for Truss Elements
# Structures).

TrussMe := module()

  # TODO: Hide module content with dummy procedures as:
  #   procname := proc(inputs); _procname(inputs); end proc:
  # so the source code is not visible by showstat and showsource commands.

  description "A Maple Library for Truss Elements Structures.";

  option package,
         load   = ModuleLoad,
         unload = ModuleUnload;

  global ground;

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
  local m_VeilingLabel;
  local m_StoredData;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Info := proc()

    description "Print module information.";

    printf(
      "+--------------------------------------------------------------------------+\n"
      "| 'TrussMe' module version 0.0 - BSD 3-Clause License - Copyright (c) 2023 |\n"
      "| Current version authors:                                                 |\n"
      "|   Matteo Larcher and Davide Stocco.                                      |\n"
      "+--------------------------------------------------------------------------+\n"
    );
    return NULL;
  end proc: # Info

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ModuleLoad := proc()

    description "Module load procedure.";

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

  export ModuleUnload := proc()

    description "Module unload procedure.";

    TrussMe:-Unprotect();
    return NULL;
  end proc: # ModuleUnload

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export TypeRegister := proc()

    description "Register module types.";

    TypeTools:-AddType('EARTH', TrussMe:-IsEarth);
    TypeTools:-AddType('FRAME', TrussMe:-IsFrame);
    TypeTools:-AddType('POINT', TrussMe:-IsPoint);
    TypeTools:-AddType('VECTOR', TrussMe:-IsVector);
    TypeTools:-AddType('BEAM', TrussMe:-IsBeam);
    TypeTools:-AddType('ROD', TrussMe:-IsRod);
    TypeTools:-AddType('RIGID_BODY', TrussMe:-IsRigidBody);
    TypeTools:-AddType('FORCE', TrussMe:-IsForce);
    TypeTools:-AddType('MOMENT', TrussMe:-IsMoment);
    TypeTools:-AddType('QFORCE', TrussMe:-IsQForce);
    TypeTools:-AddType('QMOMENT', TrussMe:-IsQMoment);
    TypeTools:-AddType('SUPPORT', TrussMe:-IsSupport);
    TypeTools:-AddType('JOINT', TrussMe:-IsJoint);
    TypeTools:-AddType('MATERIAL', TrussMe:-IsMaterial);
    TypeTools:-AddType('STRUCTURE', TrussMe:-IsStructure);
    return NULL;
  end proc: # TypeRegister

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export InitTrussMe := proc()

    description "Initialize module internal variables.";

    # Global variables
    ground := LinearAlgebra:-IdentityMatrix(4);

    # Local variables
    TrussMe:-InitLAST();
    m_gravity               := [0, 0, 0];
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
    m_VeilingLabel          := "_V";
    m_StoredData            := [];
    m_earth := table({
      "type"                 = EARTH,
      "name"                 = "earth",
      "length"               = 0,
      "frame"                = ground,
      "admissible_loads"     = [1, 1, 1, 1, 1, 1]
      }):

    return NULL;
  end proc: # InitTrussMe

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Protect := proc()

    description "Protect module variables.";

    protect(
      # Global variables
      'ground',
      # Types
      'EARTH',
      'FRAME',
      'POINT',
      'VECTOR',
      'BEAM',
      'ROD',
      'RIGID_BODY',
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

  export Unprotect := proc()

    description "Unprotect module variables.";

    unprotect(
      # Global variables
      'ground',
      # Types
      'EARTH',
      'FRAME',
      'POINT',
      'VECTOR',
      'BEAM',
      'ROD',
      'RIGID_BODY',
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

  export CheckInit := proc( $ )

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

  export InitLAST := proc(
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

  export SetLAST := proc(
    obj::LAST,
    $)

    description "Set the 'LAST' (and 'LEM') object <obj>.";

    m_LAST := obj;
    m_LEM  := m_LAST:-GetLEM(m_LAST);
    return NULL;
  end proc: # SetLAST

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetLAST := proc( $ )::LAST;

    description "Get the 'LAST' object.";

    return m_LAST;
  end proc: # GetLAST

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetLEM := proc(
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

  export `union` := proc(
    A::{list, set},
    B::{list, set},
    $)::{list, set};

    option overload;

    description "Extension of union operator to list or set objects <A> and <B>.";

    if type(A, set) and type(B, set) then
      return {op(A), op(B)};
    else
      return [op(A), op(B)];
    end if;
  end proc: # union

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetModuleOptions := proc(
    {
      VerboseMode::{integer, nothing}      := NULL,
      WarningMode::{boolean, nothing}      := NULL,
      TimeLimit::{constant, nothing}       := NULL,
      LAST_VerboseMode::{boolean, nothing} := NULL,
      LAST_WarningMode::{boolean, nothing} := NULL,
      LAST_TimeLimit::{constant, nothing}  := NULL
    }, $)

    description "Set the module options: <VerboseMode> = [0, 1, 2], <WarningMode> "
      "= [true, false], <TimeLimit> = [0, inf], <LAST_VerboseMode> = [true, false], "
      "<LAST_WarningMode> = [true, false], <LAST_TimeLimit> = [0, inf].";

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

    if (LAST_VerboseMode <> NULL) then
      if not type(LAST_VerboseMode, boolean) then
        error "invalid LAST verbose mode detected.";
      else
        m_LAST:-SetVerboseMode(m_LAST, LAST_VerboseMode);
      end if;
    end if;

    if (LAST_WarningMode <> NULL) then
      if not type(LAST_WarningMode, boolean) then
        error "invalid LAST warning mode detected.";
      else
        m_LAST:-SetWarningMode(m_LAST, LAST_WarningMode);
      end if;
    end if;

    if (LAST_TimeLimit <> NULL) then
      if (LAST_TimeLimit < 0) then
        error "invalid LAST time limit detected.";
      else
        m_LAST:-SetTimeLimit(m_LAST, LAST_TimeLimit);
      end if;
    end if;
    return NULL;
  end proc: # SetModuleOptions

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableVerboseMode::static := proc( $ )

    description "Enable the verbose mode of the module.";

    m_VerboseMode := 1;
    m_LAST:-EnableVerboseMode(m_LAST);
    return NULL;
  end proc: # EnableVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableVerboseMode::static := proc( $ )

    description "Disable the verbose mode of the module.";

    m_VerboseMode := false;
    m_LAST:-DisableVerboseMode(m_LAST);
    return NULL;
  end proc: # DisableVerboseMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export EnableWarningMode::static := proc( $ )

    description "Enable the warning mode of the module.";

    m_WarningMode := true;
    return NULL;
  end proc: # EnableWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export DisableWarningMode::static := proc( $ )

    description "Disable the warning mode of the module.";

    m_WarningMode := 0;
    return NULL;
  end proc: # DisableWarningMode

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export IsEarth := proc(
    obj::anything,
    $)::boolean;

    description "Check if the object <obj> is a EARTH object.";

    return type(obj, table) and evalb(obj["type"] = EARTH);
  end proc: # IsEarth

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export SetGravity := proc(
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

  export GetGravity := proc( $ )::Vector;

    description "Get gravity vector.";

    return m_gravity;
  end proc: # GetGravity

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Norm2 := proc(
    vec::{list, vector},
    $)::algebraic;

    description "Compute the Euclidean norm of the input vector <vec>.";

    return sqrt(add(x, x in vec^~2));
  end proc: # Norm2

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export ListPadding := proc(
    lst::{algebraic, list, Vector},
    num::integer,
    val::algebraic := 0,
    $)::{list, Vector};

    description "Pad a list or vector <lst> with value <val> to have <num> "
      "elements.";

    local out;

    if type(lst, algebraic) then
      out := [lst];
    elif type(lst, Vector) then
      out := [op(convert(lst, list))];
    else
      out := lst;
    end if;

    if (nops(out) < num) then
      out := out union [seq(val, i = 1..num-nops(out))];
    elif (nops(out) > num) then
      out := out[1..num];
    end if;

    if type(lst, Vector) then
      return convert(out, Vector);
    else
      return out;
    end if;
  end proc: # ListPadding

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Show := proc(
    tab::table,
    $)

    description "Show the content of a table <tab>.";

    print(tab = tab["type"](op(op(tab))));
    return NULL;
  end proc: # Show

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetNames := proc(
    objs::{
      list({MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH}),
      set({MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH})
    }, $)::{list({string}), set({string})};

    description "Get names of a list or set of structural objects <objs>.";

    if type(objs, set) then
      return {seq(objs[i]["name"], i = 1..nops(objs))};
    else
      return [seq(objs[i]["name"], i = 1..nops(objs))];
    end if;
  end proc: # GetNames

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetObjByName := proc(
    fld_name::string,
    objs::{
      list({MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH}),
      set({MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH})
    }, $)::{anything, MATERIAL, BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH};
    # FIXME: anything is not correct

    description "Get object which field name is <fld_name> from a list or set "
      "of objects <objs>.";

    local out, i;

    out := NULL;
    for i in objs do
      if (i["name"] = fld_name) then
        out := eval(i); # Do not remove eval: eval(table)
        break;
      end if;
    end do;
    return eval(out);
  end proc: # GetObjByName

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export GetObjsByType := proc(
    fld_type::{list(symbol), set(symbol)},
    objs::{list, set},
    $)::list;

    description "Get objects which field type is in <fld_type> from a list or "
      "set of objects <objs>.";

    local out, i;

    out := [];
    for i in objs do
      if (i::convert(fld_type, set)) then
        out := out union [eval(i)]; # Do not remove eval: eval(table)
      end if;
    end do;
    return eval(out);
  end proc: # GetObjsByType

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Simplify := proc(
    obj::anything,
    opt::anything := NULL,
    $)::anything;

    description "Simplify an algebraic expression <obj>.";

    local out, time_limit;

    time_limit := `if`(procname::indexed, op(procname), m_TimeLimit);
    try
      out := timelimit(time_limit, simplify(obj, opt));
    catch:
      WARNING("time limit of %1s exceeded for simplify operation, raw solutions "
        "is returned. The input <TimeLimit> can be modified by setting it in the "
        "TrussMe:-SetModuleOptions(...) procedure.", time_limit
      );
      out := obj;
    end try:
    return out;
  end proc: # Simplify

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Subs := proc(
    # _passed
    )::anything;

    description "Perform substitution command neglecting sub-lists and sub-sets "
      "from the substitution list.";

    map(x -> map(remove, y -> type(y, {list, set}), x), [_passed[1..-2]]);
    return subs(op(%), _passed[-1]);
  end proc: # Subs

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Diff := proc(
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
  end proc: # Diff

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export InverseFrame := proc(
    RF::FRAME,
    $)::FRAME;

    description "Inverse of the affine transformation matrix <RF>.";

    LinearAlgebra:-Transpose(RF[1..3, 1..3]);
    return <<% | -%.RF[1..3, 4]>,
            <0 | 0 | 0 | 1>>;
  end proc: # InverseFrame

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export IsFrame := proc(
    obj::anything,
    $)::boolean;

    description "Check if the object <obj> is a FRAME object.";

    return type(obj, Matrix) and
          (LinearAlgebra:-RowDimension(obj) = 4) and
          (LinearAlgebra:-ColumnDimension(obj) = 4);
  end proc: # IsFrame

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Rotate := proc(
    axis::{symbol, string},
    angle::algebraic,
    $)::FRAME;

    description "Transformation matrix corresponding to the rotation <angle> "
      "around the given <axis>.";

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

  export Translate := proc(
    x::algebraic,
    y::algebraic,
    z::algebraic,
    $)::FRAME;

    description "Affine transformation matrix corresponding to a translation "
      "<x, y, z>.";

    return <<1, 0, 0, 0>|
            <0, 1, 0, 0>|
            <0, 0, 1, 0>|
            <x, y, z, 1>>;
  end proc: # Translate

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Origin := proc(
    RF::FRAME,
    $)::Vector;

    description "Extract the origin of the reference frame <RF>.";

    return <RF[1, 4], RF[2, 4], RF[3, 4], 1>;
  end proc: # Origin

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export IsPoint := proc(
    obj::anything,
    $)::boolean;

    description "Check if the object <obj> is a POINT object.";

    return type(obj, list) and (nops(obj) = 3);
  end proc: # IsPoint

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export IsVector := proc(
    obj::anything,
    $)::boolean;

    description "Check if the object <obj> is a VECTOR object.";

    return type(obj, list) and (nops(obj) = 3);
  end proc: # IsVector

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Uvec := proc(
    axis::symbol,
    RF::FRAME := ground,
    $)::Vector;

    description "Extract the unit vector of the reference frame <RF> along the "
      "given <axis>.";

    if (axis = 'X') then
      return TrussMe:-UvecX(RF);
    elif (axis = 'Y') then
      return TrussMe:-UvecY(RF);
    elif (axis = 'Z') then
      return TrussMe:-UvecZ(RF);
    else
      error "invalid axis detected.";
    end if;
  end proc: # Uvec

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export UvecX := proc(
    RF::FRAME := ground,
    $)::Vector;

    description "Extract the x-axis unit vector of the reference frame <RF>.";

    return <RF[1, 1], RF[2, 1], RF[3, 1], 0>;
  end proc: # UvecX

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export UvecY := proc(
    RF::FRAME := ground,
    $)::Vector;

    description "Extract the y-axis unit vector of the reference frame <RF>.";

    return <RF[1, 2], RF[2, 2], RF[3, 2], 0>;
  end proc: # UvecY

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export UvecZ := proc(
    RF::FRAME := ground,
    $)::Vector;

    description "Extract the z-axis unit vector of the reference frame <RF>.";

    return <RF[1, 3], RF[2, 3], RF[3, 3], 0>;
  end proc: # UvecZ

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export Project := proc(
    x::{list, Vector},
    RF_ini::FRAME,
    RF_end::FRAME,
    $)::{list, Vector};

    description "Project the list of vector <x> (of the form [x, y, z], "
      "[x, y, z, 0], or [x, y, z, 1]) from reference frame <RF_ini> to "
      "reference frame <RF_end>.";

    local x_tmp, out, i;

    # Pad input vector with zeros
    if (nops(x) = 3) or (nops(x) = 4) then
      x_tmp := TrussMe:-ListPadding(convert(x, Vector), 4);
    else
      error "invalid input list/vector <x> detected.";
    end if;

    # Try to compare reference frames
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

  export MakeMaterial := proc({
      name::string               := "DeafultSteel",
      elastic_modulus::algebraic := 210.0E+09,
      poisson_ratio::algebraic   := 0.3,
      shear_modulus::algebraic   := elastic_modulus/(2*(1+poisson_ratio)),
      density::algebraic         := 7.4E+03
    }, $)::MATERIAL;

    description "Define a MATERIAL object with inputs: name of the material "
      "<name>, elastic modulus <elastic_modulus> (default = 210.0E9 (Pa)), "
      "Poisson's ratio <poisson_ratio> (default = 0.3 (-)), shear modulus "
      "<shear_modulus> (default = E/(2*(1+nu))), density <density> (default = "
      "7.4E3 (kg/m^3)).";

    return table({
      "type"            = MATERIAL,
      "name"            = name,
      "elastic_modulus" = elastic_modulus,
      "poisson_ratio"   = poisson_ratio,
      "shear_modulus"   = shear_modulus,
      "density"         = density
    });
  end proc: # DefineMaterial

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export IsMaterial := proc(
    obj::anything,
    $)::boolean;

    description "Check if the object <obj> is a MATERIAL object.";

    return type(obj, table) and evalb(obj["type"] = MATERIAL);
  end proc: # IsMaterial

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export MakeForce := proc(
    components::{list, Vector},
    coords::{algebraic, list(algebraic)},
    obj::{BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH},
    RF::FRAME := ground,
    $)::FORCE;

    description "Define a FORCE object with inputs: force components "
      "<components>, force application axial coordinate <coords>, target object "
      "<obj>, and optional reference frame <RF> in which the force is defined "
      "(default = ground).";

    local proj_components, admissible_components;

    if not (nops(components) = 3) then
      error "invalid input vector <components> detected.";
    end if;

    proj_components       := TrussMe:-Project(components, RF, obj["frame"]);
    admissible_components := convert(
      proj_components .~ <obj["admissible_loads"][1..3]>, list
    );
    if (proj_components <> admissible_components) and m_WarningMode then
      ["x_comp", "y_comp", "z_comp"] =~ convert(proj_components .~ <eval(
          map((x -> evalb(x = 0)), obj["admissible_loads"][1..3]),
          [true = 1, false = 0]
        )>, list);
      WARNING("Force components are not admissible for the target object. The "
        "following components will be ignored: %1", remove(x -> rhs(x) = 0, %));
    end if;

    if TrussMe:-IsSupport(obj) or TrussMe:-IsJoint(obj) then
      if (TrussMe:-ListPadding(coords, 3) <> [0, 0, 0]) then
        error "only null axial coordinate is accepted for SUPPORT and JOINT "
          "type objects.";
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

  export IsForce := proc(
    obj::anything,
    $)::boolean;

    description "Check if the object <obj> is a FORCE object.";

    return type(obj, table) and evalb(obj["type"] = FORCE);
  end proc: # IsForce

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export MakeMoment := proc(
    components::{list, Vector},
    coords::{algebraic, list(algebraic)},
    obj::{BEAM, RIGID_BODY, SUPPORT, JOINT, EARTH},
    RF::FRAME := ground,
    $)::MOMENT;

    description "Define a MOMENT object with inputs: moment components "
      "<components>, moment application axial coordinate <coords>, target object "
      "<obj>, and optional reference frame <RF> in which the moment is defined "
      "(default = ground).";

    local proj_components, admissible_components;

    if not (nops(components) = 3) then
      error "invalid input vector <components> detected.";
    end if;

    proj_components       := TrussMe:-Project(components, RF, obj["frame"]);
    admissible_components := convert(
      proj_components .~ <obj["admissible_loads"][4..6]>, list
    );
    if (proj_components <> admissible_components) and m_WarningMode then
      ["x_comp", "y_comp", "z_comp"] =~ convert(proj_components .~ <eval(
          map((x -> evalb(x = 0)), obj["admissible_loads"][4..6]),
          [true = 1, false = 0]
        )>, list);
      WARNING("Moment components are not admissible for the target object. The "
        "following components will be ignored: %1", remove(x -> rhs(x) = 0, %));
    end if;

    if (TrussMe:-IsSupport(obj) or TrussMe:-IsJoint(obj)) and
      (TrussMe:-ListPadding(coords, 3) <> [0, 0, 0]) then
      error "only null axial coordinate is accepted for SUPPORT and JOINT "
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

  export IsMoment := proc(
    obj::anything,
    $)::boolean;

    description "Check if the object <obj> is a MOMENT object.";

    return type(obj, table) and evalb(obj["type"] = MOMENT);
  end proc: # IsMoment

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export MakeQForce := proc(
    components::{procedure, list(algebraic)},
    obj::{BEAM, ROD},
    RF::FRAME := ground,
    {
      ell_min::algebraic := 0,
      ell_max::algebraic := obj["length"]
    }, $)::QFORCE;

    description "Define a QFORCE object with inputs: distributed load components "
      "<components>, target object <obj>, optional reference frame <RF> in which "
      "the load components are defined (default = ground), and optional initial "
      "<ell_min> and final <ell_max> application points (axial coordinates).";

    local proj_components, x;

    if type(components, procedure) then
      proj_components := unapply(
        TrussMe:-Project(components(x), RF, obj["frame"]), x
      );
    else
      proj_components := x -> piecewise(
        (ell_min <= x) and (x <= ell_max),
        TrussMe:-Project(components, RF, obj["frame"]),
        0
      );
    end if;

    if TrussMe:-IsRod(obj) then
      if (proj_components(x)[2] <> 0) or (proj_components(x)[3] <> 0) then
        error "only axial loads are accepted in ROD objects"
      end if;
    end if;

    return table({
      "type"       = QFORCE,
      "components" = proj_components,
      "target"     = obj["name"]
    });
  end proc: # MakeQForce

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export IsQForce := proc(
    obj::anything,
    $)::boolean;

    description "Check if the object <obj> is a QFORCE object.";

    return type(obj, table) and evalb(obj["type"] = QFORCE);
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

    description "Define a QMOMENT object with inputs: distributed torque "
      "components <components>, target object <obj>, optional reference frame "
      "<RF> in which the load components are defined (default = ground), and "
      "optional initial <ell_min> and final <ell_max> application points "
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

  export IsQMoment := proc(
    obj::anything,
    $)::boolean;

    description "Check if the object <obj> is a QMOMENT object.";

    return type(obj, table) and evalb(obj["type"] = QMOMENT);
  end proc: # IsQMoment

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export MakeSupport := proc(
    name::string,
    constrained_dof::list,
    objs::{list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})},
    coords::list,
    RF::FRAME := ground,
    {
      stiffness::{procedure,list(algebraic)} := [
        infinity, infinity, infinity,
        infinity, infinity, infinity
      ] *~ constrained_dof
    }, $)::SUPPORT;

    description "Make a SUPPORT object with inputs: support name <name>, "
      "constrained degrees of freedom <constrained_dof>, target objects <objs>, "
      "support locations <coords>, and optional reference frame <RF> in which "
      "the support is defined (default = ground). The optional input <stiffness> "
      "is a list of stiffness components (default = infinite) in the order: "
      "[ktx, kty, ktz, krx, kry, krz].";

    local S, J_tmp, i, j, sr_F_names, sr_F_values_tmp, sr_M_names, sr_M_values_tmp,
      S_stiffness, x, obj_coords;

    # Substitute -1 entries of coords with the corresponding object length
    obj_coords := [
      seq(`if`(coords[i] = -1, objs[i]["length"], coords[i]), i = 1..nops(coords))
    ];

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
      if has(map(evalb, stiffness(x)[remove(
          x -> x = 0, ([seq(i, i = 1..6)]) *~ constrained_dof
        )] <>~ 0), false) then
        error "stiffness corresponding to constrained degrees of freedom cannot be zero.";
      end if;
      # Check for zero stiffness on unconstrained dof
      if (S_stiffness(x) <> stiffness(x)) and m_WarningMode then
        WARNING("stiffness components not corresponding to constrained_dof are "
          "ignored");
      end if;
    else
      # Check for non zero stiffness on constrained dof
      if has(stiffness[remove(
          x -> x = 0, ([seq(i, i = 1..6)]) *~ constrained_dof
        )], 0) then
        error "stiffness corresponding to constrained degrees of freedom cannot "
          "be zero.";
      end if;
      # Check for zero stiffness on unconstrained dof
      S_stiffness := x -> stiffness *~ constrained_dof;
      if (S_stiffness(x) <> stiffness) and m_WarningMode then
        WARNING("stiffness components not corresponding to constrained_dof are "
          "ignored");
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
    J_tmp := TrussMe:-MakeJoint(
      name, constrained_dof, [m_earth, op(objs)], S["coordinates"], RF
    );

    S["variables"]        := J_tmp["variables"];
    S["forces"]           := J_tmp["forces"];
    S["moments"]          := J_tmp["moments"];
    S["constraint_loads"] := J_tmp["constraint_loads"];

    # Retrieve support force reactions
    sr_F_names := [FX, FY, FZ];
    for i from 1 to nops(S["forces"]) do
      if (S["forces"][i]["target"] = m_earth["name"]) then
        # Project forces in the support reference frame
        sr_F_values_tmp := TrussMe:-Project(
          S["forces"][i]["components"], ground, S["frame"]
        );
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
        sr_M_values_tmp := TrussMe:-Project(
          S["moments"][i]["components"], ground, S["frame"]
        );
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

  export IsSupport := proc(
    obj::anything,
    $)::boolean;

    description "Check if the object <obj> is a SUPPORT object.";

    return type(obj, table) and evalb(obj["type"] = SUPPORT);
  end proc: # IsSupport

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export IsCompliantSupport := proc(
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
  end proc: # IsCompliantSupport

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export CleanSupport := proc(
    obj::SUPPORT,
    $)

    description "Clean SUPPORT object <obj> internal variables.";

    obj["constraint_loads"]  := [];
    obj["support_reactions"] := [];
    obj["displacements"]     := [];
    return NULL;
  end proc: # CleanSupport

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export MakeJoint := proc(
    name::string,
    constrained_dof::list,
    objs::list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT, EARTH}),
    coords::list,
    RF::FRAME := ground,
    {
      stiffness::{procedure, list(algebraic)} := [
        infinity, infinity, infinity,
        infinity, infinity, infinity
      ] *~ constrained_dof,
      shell_objs::
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
    obj_coords := [
      seq(`if`(coords[i] = -1, objs[i]["length"], coords[i]), i = 1..nops(coords))
    ];

    for i from 1 to nops(objs) do

      if TrussMe:-IsRod(objs[i]) then
        # x coordinate of the joint location
        TrussMe:-ListPadding(eval(obj_coords[i]), 3)[1];
        # x coordinate of the joint location minus target length
        eval(eval(%^2) - eval(objs[i]["length"]^~2));
        if  %% <> 0 and %% <> 0. and
            % <> 0 and % <> 0. then
          error "JOINT objects can only be applied at extremes of ROD objects.";
        end if;
      end if;
      if TrussMe:-IsRod(objs[i]) and (constrained_dof[4..6] <> [0, 0, 0]) then
        error "ROD objects supports can only have translational constraints.";
      end if;
    end do;

    if type(stiffness, procedure) then
      J_stiffness := unapply(stiffness(x) *~ constrained_dof, x);
      # Check for non zero stiffness on constrained dof
      if has(map(evalb, stiffness(x)[remove(
          x -> x=0, ([seq(i, i = 1..6)]) *~ constrained_dof
        )] <>~ 0), false) then
        error "stiffness corresponding to constrained degrees of freedom cannot "
          "be zero.";
      end if;
      # Check for zero stiffness on unconstrained dof
      if (J_stiffness(x) <> stiffness(x)) and m_WarningMode then
        WARNING("stiffness components not corresponding to constrained_dof are "
          "ignored");
      end if;
    else
      # Check for non zero stiffness on constrained dof
      if has(stiffness[remove(
          x -> x = 0, ([seq(i, i = 1..6)]) *~ constrained_dof
        )], 0) then
        error "stiffness corresponding to constrained degrees of freedom cannot "
          "be zero.";
      end if;
      # Check for zero stiffness on unconstrained dof
      J_stiffness := x -> stiffness *~ constrained_dof;
      if (J_stiffness(x) <> stiffness) and m_WarningMode then
        WARNING("stiffness components not corresponding to constrained_dof are "
          "ignored");
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
    #      error "Joint locations are not the same on all objects.";
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
          TrussMe:-Project(jf_comp_cons, RF, objs[i]["frame"]) .~ <eval(
            map((x->evalb(x = 0)), objs[i]["admissible_loads"][1..3]),
            [true = 1, false = 0]
          )>, list);
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
        JM_||(name)||_||(objs[i]["name"]) := TrussMe:-MakeMoment(
          jm_comp_obj, obj_coords[i], objs[i], objs[i]["frame"]
          );
        # Moment on joint
        JM_||(objs[i]["name"])||_||(name) := TrussMe:-MakeMoment(
        -jm_comp_cons, 0, J, RF);
        # Use the non admissible loads to build the loads constraint
        constraint := convert(
          TrussMe:-Project(jm_comp_cons, RF, objs[i]["frame"]) .~ <eval(
            map((x->evalb(x = 0)), objs[i]["admissible_loads"][4..6]),
            [true = 1, false = 0]
          )>, list);
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

  export IsJoint := proc(
    obj::anything,
    $)::boolean;

    description "Check if the object <obj> is a JOINT object.";

    return type(obj, table) and evalb(obj["type"] = JOINT);
  end proc: # IsJoint

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export IsCompliantJoint := proc(
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
  end proc: # IsCompliantJoint

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export CleanJoint := proc(
    obj::JOINT,
    $)

    description "Clean JOINT object <obj> internal variables.";

    obj["constraint_loads"] := [];
  end proc: # CleanJoint

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export MakeRodPoints := proc(
    name::string,
    p_1::POINT,
    p_2::POINT,
    vec::{list, Vector},
    {
      area::{algebraic, procedure} := infinity,
      material::MATERIAL           := TrussMe:-MakeMaterial()
    }, $)::ROD;

    description "Create a ROD object with inputs: object name <name>, first "
      "point <p_1>, second point <p_2>, vector in XY-plane <vec>, optional "
      "section area <area> and material type <material>.";

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

    return TrussMe:-MakeRod(
      name,
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

  export MakeRod := proc(
    name::string,
    ell::algebraic,
    RF::FRAME := ground,
    {
      area::{algebraic, procedure} := infinity,
      material::MATERIAL           := TrussMe:-MakeMaterial()
    }, $)::ROD;

    description "Create a ROD object with inputs: object name <name>, reference "
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
      "name"                = name,
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

  export IsRod := proc(
    obj::anything,
    $)::boolean;

    description "Check if the object <obj> is a ROD object.";

    return type(obj, table) and evalb(obj["type"] = ROD);
  end proc: # IsRod

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export CleanRod := proc(
    obj::ROD,
    $)

    description "Clean ROD object <obj> internal variables.";

    obj["internal_actions"] := [];
    obj["displacements"]    := [];
    return NULL;
  end proc: # CleanRod

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export MakeBeamPoints := proc(
    name::string,
    p_1::POINT,
    p_2::POINT,
    vec::{list, Vector},
    {
      area::{algebraic, procedure}                   := infinity,
      timo_shear_coeff::{list(algebraic), procedure} := [5/6, 5/6],
      material::MATERIAL                             := TrussMe:-MakeMaterial(),
      I_xx::{algebraic, procedure}                   := infinity,
      I_yy::{algebraic, procedure}                   := infinity,
      I_zz::{algebraic, procedure}                   := infinity
    }, $)::BEAM;

    description "Create a BEAM object with inputs: object name <name>, first "
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
      name,
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

  export MakeBeam := proc(
    name::string,
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

    description "Create a BEAM object with inputs: object name <name>, reference "
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
      "frame_displacements" = LinearAlgebra:-IdentityMatrix(4)
    });
  end proc: # MakeBeam

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export IsBeam := proc(
    obj::anything,
    $)::boolean;

    description "Check if the object <obj> is a BEAM object.";

    return type(obj, table) and evalb(obj["type"] = BEAM);
  end proc: # IsBeam

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export CleanBeam := proc(
    obj::BEAM,
    $)

    description "Clean BEAM object <obj> internal variables.";

    obj["internal_actions"] := [];
    obj["displacements"]    := [];
    return NULL;
  end proc: # CleanBeam

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export MakeRigidBody := proc(
    name::string,
    RF::FRAME := ground,
    {
      COM::list(algebraic) := [0, 0, 0],
      mass::algebraic      := 0
    },
    $)::RIGID_BODY;

    description "Create a RIGID_BODY object with inputs: object name <name>, "
      "reference frame <RF> in which the rigid body is defined, and optional "
      "center of mass position <COM> and mass <mass>.";

    return table({
      "type"                = RIGID_BODY,
      "name"                = name,
      "frame"               = RF,
      "COM"                 = COM,
      "mass"                = mass,
      "admissible_loads"    = [1, 1, 1, 1, 1, 1],
      "frame_displacements" = LinearAlgebra:-IdentityMatrix(4)
    });
  end proc: # MakeRigidBody

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export IsRigidBody := proc(
    obj::anything,
    $)::boolean;

    description "Check if the object <obj> is a RIGID_BODY object.";

    return type(obj, table) and evalb(obj["type"] = RIGID_BODY);
  end proc: # IsRigidBody

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export MakeStructure := proc(
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
      WARNING("the structure is underconstrained with %1 unconstrained "
        "directions. Results computation may fail due to rigid body motions.",
        num_dof);
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
          WARNING("gravity load is not supported for ROD type object %1", obj);
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

  export IsStructure := proc(
    obj::anything,
    $)::boolean;

    description "Check if the object <obj> is a STRUCTURE object.";

    return type(obj, table) and evalb(obj["type"] = STRUCTURE);
  end proc: # IsStructure

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export CleanStructure := proc(
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
      elif TrussMe:-IsJoint(obj[i]) then
        obj["objects"][i] := TrussMe:-CleanJoint(i);
      end if;
    end do;
    return NULL;
  end proc: # CleanStructure

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  export CopyStructure := proc(
    struct::STRUCTURE,
    $)::STRUCTURE;

  description "Create a copy of the structure <struct> and its objects.";

    local struct_copy, obj, action;

    # Create a copy of the structure
    struct_copy := copy(struct);

    # Substitute objects in the structure with a copy
    for obj in struct_copy["objects"] do
      struct_copy["objects"] := remove(
        x -> x["name"] = obj["name"], struct_copy["objects"]
      );
      struct_copy["objects"] := struct_copy["objects"] union {copy(obj)};
    end do;

    # Substitute external actions in the structure with a copy
    for action in struct_copy["external_actions"] do
      struct_copy["external_actions"] := remove(
        x -> x["name"] = action["name"], struct_copy["external_actions"]
      );
      struct_copy["external_actions"] := struct_copy["external_actions"] union {copy(action)};
    end do;

    return struct_copy;
  end proc: # CopyStructure

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

$include "./lib/TrussMe_Plots.mpl"
$include "./lib/TrussMe_Solvers.mpl"

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

end module: # TrussMe

# That's all folks!
