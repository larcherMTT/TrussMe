
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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ObjectColor := proc(
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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotRigidBody := proc(
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

  description "Plot the RIGID_BODY object <obj> given the supports/joints "
    "list or set <joints>, the concentrated loads list or set <c_loads>, and "
    "the optional list or set of substitution <data>.";

  local p_1, p_2, i, j, idx, lines;

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
        parse("thickness") = 6
    )];
  end do;

  for j in c_loads do
    p_2 :=  subs(op(data),
      TrussMe:-Project([op(j["coordinate"]), 1], obj["frame"], ground)
    );
    lines := lines union [
      plottools:-line(
        convert(p_1[1..3], list),
        convert(p_2[1..3], list),
        parse("thickness") = 6
    )];
  end do;

  return plots:-display(
    lines,
    parse("linestyle") = solid,
    parse("color")     = TrussMe:-ObjectColor(obj),
    parse("scaling")   = constrained
  );
end proc: # PlotRigidBody

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotDeformedRigidBody := proc(
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
    data::{list(`=`), set(`=`)} := [],
    scaling::numeric            := 1
  },
  $)::function;

  description "Plot the deformed RIGID_BODY object <obj> given the "
    "supports/joints list or set <joints>, the concentrated loads list or set "
    "<c_loads>, the optional list or set of substitution <data> and  scaling "
    "factor <scaling>.";

  local p_1, p_2, rfd, js, idx, lines, j;

  lines := [];
  rfd := subs(op(data),
      obj["frame"] . ((obj["frame_displacements"] -
      LinearAlgebra:-IdentityMatrix(4)) *~ scaling +
      LinearAlgebra:-IdentityMatrix(4))
    );
  p_1 := subs(op(data), TrussMe:-Project([op(obj["COM"]), 1], rfd, ground));
  for js in joints do
    member(obj["name"], js["targets"], 'idx');
    p_2 := subs(op(data), TrussMe:-Project(
        [op(js["coordinates"][idx]), 1], rfd, ground)
      );
    lines := lines union [plottools:-line(
          convert(p_1[1..3], list), convert(p_2[1..3], list),
          parse("thickness") = 6
        )];
  end do;

  for j in c_loads do
    p_2 := subs(op(data), TrussMe:-Project(
        [op(j["coordinate"]), 1], rfd, ground)
      );
    lines := lines union [plottools:-line(
        convert(p_1[1..3], list), convert(p_2[1..3], list),
        parse("thickness") = 6
      )];
  end do;

  return plots:-display(
    lines,
    parse("linestyle") = solid,
    parse("color")     = TrussMe:-ObjectColor(obj),
    parse("scaling")   = constrained
  );
end proc: # PlotDeformedRigidBody

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotBeam := proc(
  obj::BEAM,
  {
    data::{list(`=`), set(`=`)} := []
  },
  $)::function;

  description "Plot the 'BEAM' object <obj> given a list or set of "
    "substitutions <data>.";

  local p_1, p_2;

  p_1 := subs(op(data), TrussMe:-Origin(obj["frame"]));
  p_2 := subs(op(data),
    TrussMe:-Project([obj["length"], 0, 0, 1], obj["frame"], ground)
  );

  return plots:-display(
    plottools:-line(
      convert(p_1[1..3], list),
      convert(p_2[1..3], list),
      parse("thickness") = 6
    ),
    parse("linestyle") = solid,
    parse("color")     = TrussMe:-ObjectColor(obj),
    parse("scaling")   = constrained
  );
end proc: # PlotBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotDeformedBeam := proc(
  obj::BEAM,
  {
    data::{list(`=`), set(`=`)} := [],
    scaling::numeric            := 1
  },
  $)::function;

  description "Plot the deformed 'BEAM' object <obj> given an optional list or "
    "set of substitution <data> and a scaling factor <scaling>.";

  local sc, rfd;

  rfd := subs(op(data), obj["frame"] . ((obj["frame_displacements"] -
    LinearAlgebra:-IdentityMatrix(4)) *~ scaling +
    LinearAlgebra:-IdentityMatrix(4)));

  sc := subs(op(data), TrussMe:-Project(
    subs(obj["displacements"](x),
      [x, 0, 0, 0] +~ [ux(x)*~scaling, uy(x)*~scaling, uz(x)*~scaling, 1]
      ), rfd, ground)[1..3]);

  return plots:-display(
    plots:-spacecurve(
      sc, x = subs(op(data), 0..obj["length"]),
      parse("thickness") = 6
    ),
    parse("linestyle") = solid,
    parse("color")     = TrussMe:-ObjectColor(obj),
    parse("scaling")   = constrained
  );
end proc: # PlotDeformedBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotRod := proc(
  obj::ROD,
  {
    data::{list(`=`), set(`=`)} := []
  },
  $)::function;

  description "Plot the ROD object <obj> given an optional list or set of "
    "substitution data <data>.";

  local p_1, p_2;

  p_1 := subs(op(data), TrussMe:-Origin(obj["frame"]));
  p_2 := subs(op(data),
    TrussMe:-Project([obj["length"], 0, 0, 1], obj["frame"], ground)
  );

  return plots:-display(
    plottools:-line(
      convert(p_1[1..3], list),
      convert(p_2[1..3], list),
      parse("thickness") = 4
    ),
    parse("linestyle") = solid,
    parse("color")     = TrussMe:-ObjectColor(obj),
    parse("scaling")   = constrained
  );
end proc: # PlotRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotDeformedRod := proc(
  obj::ROD,
  {
    data::{list(`=`), set(`=`)} := [],
    scaling::numeric            := 1
  },
  $)::function;

  description "Plot the ROD object <obj> given the list of <targets>, an "
    "optional list or set of substitution data <data> and scaling factor "
    "<scaling>.";

  local p_1, p_2, rfd;

  rfd := obj["frame"] . ((obj["frame_displacements"] -
    LinearAlgebra:-IdentityMatrix(4)) *~ scaling +
    LinearAlgebra:-IdentityMatrix(4));

  p_1 := subs(op(data), TrussMe:-Origin(rfd));
  p_2 := subs(op(data), TrussMe:-Project(
      [obj["length"] + subs(obj["displacements"](obj["length"]
    ), ux(obj["length"]) *~ scaling), 0, 0, 1], rfd, ground));

  return plots:-display(
    plottools:-line(
      convert(p_1[1..3], list),
      convert(p_2[1..3], list),
      parse("thickness") = 4
    ),
    parse("linestyle") = solid,
    parse("color")     = TrussMe:-ObjectColor(obj),
    parse("scaling")   = constrained
  );
end proc: # PlotDeformedRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotJoint := proc(
  obj::JOINT,
  targets::{
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  {
    data::{list(`=`), set(`=`)} := []
  },
  $)::function;

  description "Plot the JOINT object <obj> given the list of <targets> and "
    "an optional list or set of substitution data <data>.";

  local O;

  O := subs(op(data), TrussMe:-Origin(
    TrussMe:-GetObjByName(obj["targets"][1], targets)["frame"].
    TrussMe:-Translate(op(TrussMe:-ListPadding(obj["coordinates"][1],3)))
    ));

  return plots:-display(
    plottools:-point(
      convert(O[1..3], list),
      parse("symbol")     = solidsphere,
      parse("symbolsize") = 20
    ),
    parse("linestyle") = solid,
    parse("color")     = TrussMe:-ObjectColor(obj),
    parse("scaling")   = constrained
  );
end proc: # PlotJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotDeformedJoint := proc(
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

  description "Plot the JOINT object <obj> given the list of <targets>, "
    "an optional list or set of substitution data <data> and scaling factor "
    "<scaling>.";

  local O, rfd;

  # TODO: add compliant joint deformation

  rfd := (obj["frame_displacements"] - LinearAlgebra:-IdentityMatrix(4)) *~
    scaling + LinearAlgebra:-IdentityMatrix(4);

  # FIXME: joint target may be a joint itself
  O := subs(op(data), TrussMe:-Origin(
    TrussMe:-GetObjByName(obj["targets"][1], targets)["frame"] .
    TrussMe:-Translate(op(TrussMe:-ListPadding(obj["coordinates"][1],3))))[1..3] +~
    TrussMe:-Project(TrussMe:-Origin(rfd)[1..3], obj["frame"], ground)
    );

  return plots:-display(
    plottools:-point(
      convert(O, list),
      symbol     = solidsphere,
      symbolsize = 20
    ),
    parse("linestyle") = solid,
    parse("color")     = TrussMe:-ObjectColor(obj),
    parse("scaling")   = constrained
  );
end proc: # PlotDeformedJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotSupport := proc(
  obj::SUPPORT,
  targets::{
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  {
    data::{list(`=`), set(`=`)} := []
  },
  $)::function;

  description "Plot the SUPPORT object <obj> given the list of <targets>, "
    "an optional list or set of substitution data <data> and scaling factor "
    "<scaling>.";

  local O;

  if (nops(obj["targets"]) > 1) then
    O := subs(op(data), TrussMe:-Origin(
      TrussMe:-GetObjByName(obj["targets"][2], targets)["frame"].
      TrussMe:-Translate(op(TrussMe:-ListPadding(obj["coordinates"][2], 3)))
      ));
  else
    O := subs(op(data), TrussMe:-Origin(
      m_earth["frame"].
      TrussMe:-Translate(op(TrussMe:-ListPadding(obj["coordinates"][1], 3)))
      ));
  end if;

  return plots:-display(
    plottools:-point(
      convert(O[1..3], list),
      parse("symbol")     = solidbox,
      parse("symbolsize") = 20
    ),
    parse("linestyle") = solid,
    parse("color")     = TrussMe:-ObjectColor(obj),
    parse("scaling")   = constrained
  );
end proc: # PlotSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotDeformedSupport := proc(
  obj::SUPPORT,
  targets::{
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  {
    data::{list(`=`), set(`=`)} := [],
    scaling::numeric            := 1
  },
  $)::function;

  description "Plot the deformed SUPPORT object <obj> given the list of "
    "<targets>, an optional list or set of substitution data <data> and scaling "
    "factor <scaling>.";

  local O;

  # TODO: add compliant support deformation

  if (nops(obj["targets"]) > 1) then
    O := subs(op(data), TrussMe:-Origin(
      TrussMe:-GetObjByName(obj["targets"][2], targets)["frame"].
      TrussMe:-Translate(op(TrussMe:-ListPadding(obj["coordinates"][2], 3)))
      ));
  else
    O := subs(op(data), TrussMe:-Origin(
      m_earth["frame"].
      TrussMe:-Translate(op(TrussMe:-ListPadding(obj["coordinates"][1], 3)))
      ));
  end if;

  return plots:-display(
    plottools:-point(
      convert(O[1..3], list),
      parse("symbol")     = solidbox,
      parse("symbolsize") = 20
    ),
    parse("linestyle") = solid,
    parse("color")     = TrussMe:-ObjectColor(obj),
    parse("scaling")   = constrained
  );
end proc: # PlotDeformedSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotStructure := proc(
  str::STRUCTURE,
  {
    data::{list(`=`), set(`=`)} := []
  },
  $)::{function, list(function)};

  description "Plot the STRUCTURE object <str> given an optional list or set "
    "of substitution data <data>.";

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
      map(TrussMe:-GetObjByName, obj["targets"], str["objects"]);
      disp := disp union [TrussMe:-PlotSupport(obj, %, parse("data") = data)];
    elif TrussMe:-IsJoint(obj) then
      map(TrussMe:-GetObjByName, obj["targets"], str["objects"]);
      disp := disp union [TrussMe:-PlotJoint(obj, %, parse("data") = data)];
    elif TrussMe:-IsRigidBody(obj) then
      TrussMe:-GetObjsByType([JOINT, SUPPORT], str["objects"]);
      rb_joints := remove(x -> (not member(obj["name"], x["targets"])), %);
      TrussMe:-GetObjsByType([FORCE, MOMENT], str["external_actions"]);
      rb_loads := remove(x -> obj["name"] <> x["target"], %);
      disp := disp union [
        TrussMe:-PlotRigidBody(obj, rb_joints, rb_loads, parse("data") = data)
      ];
    end if;
  end do;

  return plots:-display(
    disp,
    parse("axes")    = boxed,
    parse("scaling") = constrained,
    parse("labels")  = ['x', 'y', 'z']
  );
end proc: # PlotStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export PlotDeformedStructure := proc(
  str::STRUCTURE,
  {
    data::{list(`=`), set(`=`)} := [],
    scaling::numeric            := 1
  },
  $)::{function, list(function)};

  description "Plot the deformed 'STRUCTURE' object <str> given a optional list "
    "or set of substitution data <data> and scaling factor <scaling>.";

  local disp, rb_joints, rb_loads, obj;

  # Check if displacements and frame displacements are solved
  if not str["displacements_solved"] or not str["frame_displacements_solved"] then
    error "displacements and frame displacements must be solved before plotting "
      "the deformed structure.";
  end if;

  disp := []:
  for obj in str["objects"] do
    if TrussMe:-IsBeam(obj) then
      disp := disp union [TrussMe:-PlotDeformedBeam(
          obj,
          parse("data")    = data,
          parse("scaling") = scaling
        )];
    elif TrussMe:-IsRod(obj) then
      disp := disp union [TrussMe:-PlotDeformedRod(
          obj,
          parse("data")    = data,
          parse("scaling") = scaling
        )];
    elif TrussMe:-IsSupport(obj) then
      disp := disp union [TrussMe:-PlotDeformedSupport(
          obj, map(TrussMe:-GetObjByName, obj["targets"], str["objects"]),
          parse("data")    = data,
          parse("scaling") = scaling
        )];
    elif TrussMe:-IsJoint(obj) then
      disp := disp union [TrussMe:-PlotDeformedJoint(
          obj, map(TrussMe:-GetObjByName, obj["targets"], str["objects"]),
          parse("data")    = data,
          parse("scaling") = scaling
        )];
    elif TrussMe:-IsRigidBody(obj) then
      TrussMe:-GetObjsByType([JOINT, SUPPORT], str["objects"]);
      rb_joints := remove(x -> (not member(obj["name"], x["targets"])), %);
      TrussMe:-GetObjsByType([FORCE, MOMENT], str["external_actions"]);
      rb_loads := remove(x -> obj["name"] <> x["target"], %);
      disp := disp union [TrussMe:-PlotDeformedRigidBody(
          obj, rb_joints, rb_loads,
          parse("data")    = data,
          parse("scaling") = scaling
        )];
    end if;
  end do;

  return plots:-display(
    disp,
    parse("axes")    = boxed,
    parse("scaling") = constrained,
    parse("labels")  = ['x', 'y', 'z']
  );
end proc: # PlotDeformedStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsInsideJoint := proc(
  obj::JOINT,
  pnt::POINT,
  tol::numeric := 1e-4,
  $)::boolean;

  description "Check if the point <pnt> is inside the JOINT <obj> within a "
    "optianl tolerance value <tol>.";

  local O;

  if (nops(pnt) <> 3) then
    error "the input point must be a list of 3 elements.";
  end if;

  if (nops(obj["targets"]) > 1) then
    O := TrussMe:-Origin(
      TrussMe:-GetObjByName(obj["targets"][1], targets)["frame"] .
      TrussMe:-Translate(obj["coordinates"][1], 0, 0)
      );
  elif m_WarningMode then
    WARNING("the support has no targets.");
  end if;

  return evalb(norm(pnt - O) <= tol);
end proc: # IsInsideJoint

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsInsideSupport := proc(
  obj::SUPPORT,
  pnt::POINT,
  tol::numeric := 1e-4,
  $)::boolean;

  description "Check if the point <p> is inside the SUPPORT <obj> within a "
    "optianol tolerance <tol> value.";

  local O;

  if (nops(pnt) <> 3) then
    error "the input point must be a list of 3 elements.";
  end if;

  if (nops(obj["targets"]) > 1) then
    O := TrussMe:-Origin(
      TrussMe:-GetObjByName(obj["targets"][2], targets)["frame"] .
      TrussMe:-Translate(obj["coordinates"][2], 0, 0)
      );
  elif m_WarningMode then
    WARNING("the support has no targets.");
  end if;

  return evalb(norm(pnt - O) <= tol);
end proc: # IsInsideSupport

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsInsideRod := proc(
  obj::ROD,
  pnt::POINT,
  $)::boolean;

  description "Check if the point <pnt> is inside the ROD <obj>.";

  local O, V, W;

  if (nops(pnt) <> 3) then
    error "the input point must be a list of 3 elements.";
  end if;

  O := TrussMe:-Origin(obj["frame"]);
  V := obj["frame"].TrussMe:-Translate(obj["length"], 0, 0) - O;
  W := pnt - O;
  return evalb(dot(W, V) >= 0) and (dot(W, V) <= dot(V, V));
end proc: # IsInsideRod

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsInsideBeam := proc(
  obj::BEAM,
  pnt::POINT,
  $)::boolean;

  description "Check if the point <pnt> is inside the 'BEAM' <obj>.";

  local O, V, W, out;

  if (nops(pnt) <> 3) then
    error "the input point must be a list of 3 element.";
  end if;

  O := TrussMe:-Origin(obj["frame"]);
  V := obj["frame"].TrussMe:-Translate(obj["length"], 0, 0) - O;
  W := pnt - O;
  return evalb((dot(W, V) >= 0) and (dot(W, V) <= dot(V, V)));
end proc: # IsInsideBeam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export IsInsideStructure := proc(
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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
