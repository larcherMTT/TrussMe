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

export ComputeDOF := proc(
  objs::{
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT}),
    set( {BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  }, $)::integer, function;

  description "Compute the degree of freedom of the input structure objects "
    "<objs>.";

  local dof, objs_tmp, obj, i, j, k, vertex, colors, G;

  objs_tmp := objs union [m_earth];

  # Built connections graph
  vertex := [];
  colors := [];
  if (m_VerboseMode > 0) then
    printf("TrussMe:-ComputeDOF(...): checking structure connections...\n");
  end if;
  for i from 1 to nops(objs_tmp) do
    vertex := vertex union [objs_tmp[i]["name"]];
    colors := colors union [TrussMe:-ObjectColor(objs_tmp[i])];
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
    printf("TrussMe:-ComputeDOF(...): checking structure connections... DONE\n");
  end if;

  # Check graph connections
  if GraphTheory:-IsConnected(G) then
    if (m_VerboseMode > 1) then
      print(GraphTheory:-DrawGraph(G), layout = tree);
    end if;
  else
    print(GraphTheory:-DrawGraph(G), layout = tree);
    WARNING("unconnected elements detected in the structure.");
  end if;

  if (m_VerboseMode > 0) then
    printf("TrussMe:-ComputeDOF(...): computing DOF...\n");
  end if;

  # Compute dofs
  dof := 0;
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
    printf("TrussMe:-ComputeDOF(...): computing DOF... DONE (DOF = %d)\n", dof);
  end if;

  return dof, G;
end proc: # ComputeDOF

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DrawStructureGraph := proc(
  obj::STRUCTURE,
  $)::function;

  description "Draw the connections graph of the STRUCTURE object <obj>.";

  return plots:-display(
    GraphTheory:-DrawGraph(obj["connections_graph"], layout = tree),
    title = "Structure connections graph"
  );
end proc: # DrawStructureGraph

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export DrawStructureSparseMatrix := proc(
  obj::STRUCTURE,
  {
    gauss_elimin::boolean := false
  }, $)::function;

  description "Draw the sparse matrix for the equation system of STRUCTURE "
    "object <obj> and optionally apply Gaussian elimination to the matrix with "
    "the option <gauss_elimin>.";

  local A, B;

  A, B := LinearAlgebra:-GenerateMatrix(obj["equations"], obj["variables"]);
  if gauss_elimin then
    A := LinearAlgebra:-GaussianElimination(A);
  end if;

  return plots:-display(
    plots:-sparsematrixplot(A, matrixview),
    title = "Structure sparse matrix"
  );
end proc: # DrawStructureSparseMatrix

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export NewtonEuler := proc(
  exts::{
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set({FORCE, MOMENT, QFORCE, QMOMENT})
  },
  obj::{BEAM, ROD, RIGID_BODY, SUPPORT, JOINT},
  {
    pole::POINT          := [0, 0, 0],
    upper_lim::algebraic := obj["length"]
  }, $)::list;

  description "Compute the Newton-Euler static equilibrium equations given a set "
    "of external actions <exts>, and object to compute the equilibrium <obj>, "
    "the axial coordinate of the pole <pole>, and an optional upper limit of the "
    "integration <upper_lim>.";

  local eq_T, eq_R, ext, arm;

  eq_T := [0, 0, 0];
  for ext in exts do
    if ext["target"] = obj["name"] then
      if TrussMe:-IsForce(ext) then
        eq_T := eq_T + ext["components"];
      elif TrussMe:-IsQForce(ext) then
        eq_T := eq_T + convert(map(
          integrate, ext["components"](x), x = 0..upper_lim
        ), list);
      end if;
    elif m_WarningMode then
      WARNING("TrussMe:-NewtonEuler(...): %1 is not applied to %2.", ext, obj);
    end if;
  end do;

  eq_R := [0, 0, 0];
  for ext in exts do
    if (ext["target"] = obj["name"]) then
      if TrussMe:-IsMoment(ext) then
        eq_R := eq_R + ext["components"];
      elif TrussMe:-IsForce(ext) then
        arm := convert(pole - ext["coordinate"], Vector);
        eq_R := eq_R + convert(LinearAlgebra:-CrossProduct(
          convert(ext["components"], Vector), arm
        ), list);
      elif TrussMe:-IsQForce(ext) then
        arm := convert(pole - [x, 0, 0], Vector);
        eq_R := eq_R + map(integrate, convert(LinearAlgebra:-CrossProduct(
          convert(ext["components"](x), Vector), arm
        ), list), x = 0..upper_lim);
      elif TrussMe:-IsQMoment(ext) then
        eq_R := eq_R + map(integrate, ext["components"](x), x = 0..upper_lim);
      end if;
    elif m_WarningMode then
      WARNING("TrussMe:-NewtonEuler(...): %1 is not applied to %2.", ext, obj);
    end if;
  end do;

  return eq_T union eq_R;
end proc: # NewtonEuler

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export SolveStructure := proc(
  struct::STRUCTURE,
  {
    compute_internal_actions::boolean    := false,
    compute_displacements::boolean       := false,
    compute_potential_energy::boolean    := false,
    compute_frame_displacements::boolean := false,
    timoshenko_beam::boolean             := false,
    implicit::boolean                    := false,
    unveil_results::boolean              := true,
    dummy_vars::{list, set}              := []
  }, $)::STRUCTURE;

  description "Solve the static equilibrium of a structure with inputs: "
    "structure <struct>, optional compute internal action enabling flag "
    "<compute_internal_actions>, optional compute displacement enabling flag "
    "<compute_displacements>, optional Timoshenko beam flag <timoshenko_beam>, "
    "optional implicit solution flag <implicit>, and optional unveil results "
    "flag <unveil_results>.";

  local g_load, S_obj, S_rigid, S_ext, S_support, S_joint, S_con_forces, vars,
    sol, obj, x, str_eq, str_vars, P_energy, veiling_idx, veiling_label, veils,
    i, uveils;

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

  # Add dummy_vars_subs to m_StoredData to avoid LAST pivoting on dummy variables
  dummy_vars_subs := dummy_vars =~ [seq(0, i = 1..nops(dummy_vars))];
  m_StoredData := m_StoredData union dummy_vars_subs;

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
      printf("TrussMe:-SolveStructure(...): solving the isostatic structure...\n");
    end if;
    sol, veils, str_eq, str_vars := TrussMe:-IsostaticSolver(
      S_obj union S_rigid union S_joint union S_support,
      S_ext union S_con_forces,
      vars,
      parse("implicit") = implicit
      );
    # Update Structure equations and variables
    struct["equations"] := str_eq;
    struct["variables"] := str_vars;
    if (m_VerboseMode > 0) then
      printf("TrussMe:-SolveStructure(...): solving the isostatic structure... DONE\n");
      printf("TrussMe:-SolveStructure(...): updating support reactions fields...\n");
    end if;
    # Update support reactions properties
    for obj in S_support do
      obj["support_reactions"] := subs(sol, obj["generic_support_reactions"]);
    end do;
    if (m_VerboseMode > 0) then
      printf("TrussMe:-SolveStructure(...): updating support reactions fields... DONE\n");
    end if;

  # Solve hyperstatic structure
  elif (struct["dof"] < 0) then

    if (nops(struct["hyperstatic_variables"]) <> -struct["dof"]) then
      error "TrussMe:-SolveStructure(...): mismatch in the structure degrees of freedom, check the hyper"
        "static variables of the structure and update the structure object";
    end if;
    if (m_VerboseMode > 0) then
      printf("TrussMe:-SolveStructure(...): solving the hyperstatic structure...\n");
    end if;
    sol, veils, str_eq, str_vars, P_energy := TrussMe:-HyperstaticSolver(
      S_obj union S_rigid union S_joint union S_support,
      S_ext union S_con_forces,
      vars,
      struct["hyperstatic_variables"],
      struct["hyperstatic_displacements"],
      parse("timoshenko_beam") = timoshenko_beam,
      parse("implicit")        = implicit
      );
    if (m_VerboseMode > 0) then
      printf("TrussMe:-SolveStructure(...): solving the hyperstatic structure... DONE\n");
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
      printf("TrussMe:-SolveStructure(...): updating support reactions fields...\n");
    end if;
    for obj in S_support do
    obj["support_reactions"] := subs(sol, obj["generic_support_reactions"]);
    end do;
    if (m_VerboseMode > 0) then
      printf("TrussMe:-SolveStructure(...): updating support reactions fields... "
        "DONE\n");
    end if;
  end if;

  # Add veils
  struct["veils"] := struct["veils"] union veils;

  # Set support reactions solved flag
  struct["support_reactions_solved"] := true;

  # Compute internal actions
  if (compute_internal_actions or compute_displacements or compute_potential_energy) and
     not struct["internal_actions_solved"] then
    TrussMe:-ComputeInternalActions(
      S_obj, S_ext union S_con_forces, sol
    );
    # Set internal actions computed flag
    struct["internal_actions_solved"] := true;
  end if;

  # Compute potential energy
  if compute_potential_energy and (not struct["potential_energy_solved"] or nops(dummy_vars)>0) then
    if implicit then
      error "TrussMe:-SolveStructure(...): potential energy cannot be computed "
        "in implicit mode";
    end if;
    if (m_VerboseMode > 0) then
      printf("TrussMe:-SolveStructure(...): computing potential energy...\n");
    end if;
    P_energy, uveils := TrussMe:-ComputePotentialEnergy(
      S_obj union S_support union S_joint, sol,
      parse("timoshenko_beam") = timoshenko_beam,
      parse("dummy_vars")      = dummy_vars,
      parse("veils")           = veils);
    if (m_VerboseMode > 0) then
      printf("TrussMe:-SolveStructure(...): computing potential energy... "
        "DONE\n");
    end if;
    # Update structure energy
    struct["potential_energy"] := P_energy;
    # Update structure veils
    struct["veils"] := struct["veils"] union uveils;
    # Set potential energy computed flag
    struct["potential_energy_solved"] := true;
  end if;

  # Compute displacements
  if compute_displacements and not struct["displacements_solved"] then
    TrussMe:-ComputeDisplacements(
      S_obj union S_joint union S_support, S_ext union S_con_forces, sol,
      parse("timoshenko_beam") = timoshenko_beam
    );
    # Set displacements computed flag
    struct["displacements_solved"] := true;
  end if;

  if compute_frame_displacements and not struct["frame_displacements_solved"] then
    if (m_VerboseMode > 0) then
      printf("TrussMe:-SolveStructure(...): computing frame displacements...\n");
    end if;
    veils := TrussMe:-ComputeObjectFrameDisplacements(
      struct,
      parse("timoshenko_beam") = timoshenko_beam,
      parse("unveil_results")  = unveil_results
    );
    if (m_VerboseMode > 0) then
      printf("TrussMe:-SolveStructure(...): computing frame displacements... "
        "DONE\n");
    end if;
    # Set frame displacements computed flag
    struct["frame_displacements_solved"] := true;

    # Add veils
    struct["veils"] := struct["veils"] union veils;
  end if;

  return struct;
end proc: # SolveStructure

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export HyperstaticSolver := proc(
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
    S_objs, E_objs, sol, iso_veils, hyper_veils, veils;

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
  iso_sol, iso_veils, iso_eq, iso_vars := TrussMe:-IsostaticSolver(
    objs, exts, iso_vars, parse("implicit") = implicit
  );

  # Compute internal actions
  TrussMe:-ComputeInternalActions(S_objs, exts, iso_sol);

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
  P_energy := TrussMe:-Simplify(P_energy);

  #print("P_energy: ", P_energy); #REMOVE
  #print("hyper_vars: ", hyper_vars); #REMOVE
  #print("veils: ", veils); #REMOVE
  #print("hyper_disp: ", hyper_disp); #REMOVE

  # Compute the hyperstatic equation
  hyper_eq := TrussMe:-Diff~(
    P_energy, hyper_vars, parse("veils") = iso_veils
  ) =~ hyper_disp;

  # Substitute Float(undefined) with 0 in the derivative of the potential energy
  # (this comes in case of non derivable piecewise functions in the compliant
  # joint stiffness)
  # FIXME: this is a temporary fix
  hyper_eq := TrussMe:-Simplify(eval(hyper_eq, Float(undefined) = 0));

  # Check for implicit solution flag
  if implicit then
    hyper_sol := iso_sol;
  else
    if (m_VerboseMode > 0) then
      printf("TrussMe:-HyperstaticSolver(...): solving the hyperstatic variables...\n");
    end if;
    # Solve hyperstatic equations
    #print("hyper_eq: ", hyper_eq); #REMOVE
    #print("hyper_vars: ", hyper_vars); #REMOVE
    #hyper_sol := op(solve(hyper_eq, hyper_vars)) assuming real; # Maple solver #REMOVE
    # Solve structure equations (LinearSolver)
    if m_LinearSolver = "LAST" then
      hyper_sol, hyper_veils := TrussMe:-LinearSolver(hyper_eq, hyper_vars);
    elif m_LinearSolver = "Maple" then
      hyper_sol := op(solve(hyper_eq, hyper_vars)) assuming real;
      hyper_veils := [];
    else
      error "TrussMe:-HyperstaticSolver(...): invalid linear solver.";
    end if;

    if (hyper_sol = NULL) then
      error "TrussMe:-HyperstaticSolver(...): hyperstatic solution not found.";
    end if;

    # Substitute hyper_sol in P_energy
    P_energy := subs(hyper_sol, P_energy);

    # Update m_StoredData
    m_StoredData := subs(sol, m_StoredData);

    if (m_VerboseMode > 0) then
      printf("TrussMe:-HyperstaticSolver(...): solving the hyperstatic variables... "
        "DONE\n");
    end if;
    sol := hyper_sol union subs(hyper_sol, iso_sol);
    veils := subs(hyper_sol, iso_veils) union hyper_veils;
  end if;

  if (_nresults = 5) then
    return sol, veils, iso_eq union hyper_eq, iso_vars union hyper_vars, P_energy;
  else
    return sol, veils;
  end
end proc: # HyperstaticSolver

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputePotentialEnergy := proc(
  objs::{
    list({BEAM, ROD, SUPPORT, JOINT}),
    set({BEAM, ROD, SUPPORT, JOINT})
  },
  sol::{list, set} := [],
  {
    timoshenko_beam::boolean := false,
    dummy_vars::{list, set}  := [],
    veils::{list, set}       := []
  }, $)

  description "Compute the internal potential energy of the structure given the "
    "objects <objs> and optional Timoshenko beam flag <timoshenko_beam>.";

  local dummy_vars_subs, obj, P, x, f, FJX, FJY, FJZ, MJX, MJY, MJZ, i, uveils_subs, uveils;

  dummy_vars_subs := dummy_vars =~ [seq(0, i = 1..nops(dummy_vars))];
  # Compute undummy veils
  print("here1");
  uveils_subs := lhs~(veils) =~ cat~(lhs~(veils),__ud);
  print("here1");
  subs(dummy_vars_subs, veils);
  print("here1");
  #lhs~(%) =~ TrussMe:-Simplify(eval['recurse'](rhs~(%), %));
  v1 := %;
  v_prec := [];
  while evalb(v1 <> v_prec) do
    v_prec := v1;
    v1 := lhs~(v1) =~ eval(rhs~(v1),v1);
  end do:
  v1;
  print("here1");
  uveils := subs(uveils_subs, %);
  print("here1");
  uveils_zero, uveils := selectremove(x -> rhs(x) = 0, uveils);
print("here1");
  # Update StoredData
  subs(op(m_StoredData), uveils);
  print("here1");
  m_StoredData := m_StoredData union (lhs~(%) =~ subs(op(ListTools:-Reverse(%)), rhs~(%)));
print("here1");
  P := 0;
  for obj in objs do
    if TrussMe:-IsBeam(obj) or TrussMe:-IsRod(obj) then
      # Normal action N contribution
      if (member(N, map(lhs, obj["internal_actions"]))) and
          (subs(obj["internal_actions"](x), N(x)) <> 0) then
        subs(obj["internal_actions"](x), N(x));
        subs(uveils_subs, uveils_zero, dummy_vars_subs, %);
        # (%% - %) -> internal action of the structure due to dummy loads only
        # %        -> internal action of the structure without dummy loads
        # %%       -> internal action of the structure with all loads
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
          subs(uveils_subs, uveils_zero, dummy_vars_subs, %);
          # (%% - %) -> internal action of the structure due to dummy loads only
          # %        -> internal action of the structure without dummy loads
          # %%       -> internal action of the structure with all loads
          P := P + integrate(
              eval(`if`(nops(dummy_vars) > 0, 2 * (%% - %) * %, %%^2))/
                (2*obj["timo_shear_coeff"](x)[1]*obj["material"]["shear_modulus"]*obj["area"](x)),
              x = 0..obj["length"]);
        end if;
        # Shear action Tz contribution
        if (member(Tz, map(lhs, obj["internal_actions"]))) and
            (subs(obj["internal_actions"](x), Tz(x)) <> 0) then
          subs(obj["internal_actions"](x), Tz(x));
          subs(uveils_subs, uveils_zero, dummy_vars_subs, %);
          # (%% - %) -> internal action of the structure due to dummy loads only
          # %        -> internal action of the structure without dummy loads
          # %%       -> internal action of the structure with all loads
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
        subs(uveils_subs, uveils_zero, dummy_vars_subs, %);
        # (%% - %) -> internal action of the structure due to dummy loads only
        # %        -> internal action of the structure without dummy loads
        # %%       -> internal action of the structure with all loads
        P := P + integrate(
            eval(`if`(nops(dummy_vars) > 0, 2 * (%% - %) * %, %%^2))/
              (2*obj["material"]["shear_modulus"]*obj["inertias"][1](x)),
            x = 0..obj["length"]);
          end if;
      # Bending moment action My contribution
      if (member(My, map(lhs, obj["internal_actions"]))) and
          (subs(obj["internal_actions"](x), My(x)) <> 0) then
        subs(obj["internal_actions"](x), My(x));
        subs(uveils_subs, uveils_zero, dummy_vars_subs, %);
        # (%% - %) -> internal action of the structure due to dummy loads only
        # %        -> internal action of the structure without dummy loads
        # %%       -> internal action of the structure with all loads
        P := P + integrate(
            eval(`if`(nops(dummy_vars) > 0, 2 * (%% - %) * %, %%^2))/
              (2*obj["material"]["elastic_modulus"]*obj["inertias"][2](x)),
            x = 0..obj["length"]);
      end if;
      # Bending moment action Mz contribution
      if (member(Mz, map(lhs, obj["internal_actions"]))) and
          (subs(obj["internal_actions"](x), Mz(x)) <> 0) then
        subs(obj["internal_actions"](x), Mz(x));
        subs(uveils_subs, uveils_zero, dummy_vars_subs, %);
        # (%% - %) -> internal action of the structure due to dummy loads only
        # %        -> internal action of the structure without dummy loads
        # %%       -> internal action of the structure with all loads
        P := P + integrate(
            eval(`if`(nops(dummy_vars) > 0, 2 * (%% - %) * %, %%^2))/
              (2*obj["material"]["elastic_modulus"]*obj["inertias"][3](x)),
            x = 0..obj["length"]);
      end if;
    elif TrussMe:-IsCompliantSupport(obj) then
      # Support reaction Fx contribution
      if (subs(obj["generic_support_reactions"], FX) <> 0) and
          (obj["stiffness"](x)[1] <> infinity) and
          (obj["constrained_dof"][1] <> 0) then
        P := P + subs(obj["generic_support_reactions"], sol,
            TrussMe:-ComputeSpringEnergy(-FX, (x -> obj["stiffness"](x)[1]))
          );
      end if;
      # Support reaction Fy contribution
      if (subs(obj["generic_support_reactions"], FY) <> 0) and
          (obj["stiffness"](x)[2] <> infinity) and
          (obj["constrained_dof"][2] <> 0) then
        P := P + subs(obj["generic_support_reactions"], sol,
            TrussMe:-ComputeSpringEnergy(-FY, (x -> obj["stiffness"](x)[2]))
          );
      end if;
      # Support reaction Fz contribution
      if (subs(obj["generic_support_reactions"], FZ) <> 0) and
          (obj["stiffness"](x)[3] <> infinity) and
          (obj["constrained_dof"][3] <> 0) then
        P := P + subs(obj["generic_support_reactions"], sol,
            TrussMe:-ComputeSpringEnergy(-FZ, (x -> obj["stiffness"](x)[3]))
          );
      end if;
      # Support reaction Mx contribution
      if (subs(obj["generic_support_reactions"], MX) <> 0) and
          (obj["stiffness"](x)[4] <> infinity) and
          (obj["constrained_dof"][4] <> 0) then
        P := P + subs(obj["generic_support_reactions"], sol,
            TrussMe:-ComputeSpringEnergy(-MX, (x -> obj["stiffness"](x)[4]))
          );
      end if;
      # Support reaction My contribution
      if (subs(obj["generic_support_reactions"], MY) <> 0) and
          (obj["stiffness"](x)[5] <> infinity) and
          (obj["constrained_dof"][5] <> 0) then
        P := P + subs(obj["generic_support_reactions"], sol,
            TrussMe:-ComputeSpringEnergy(-MY, (x -> obj["stiffness"](x)[5]))
          );
      end if;
      # Support reaction Mz contribution
      if (subs(obj["generic_support_reactions"], MZ) <> 0) and
          (obj["stiffness"](x)[6] <> infinity) and
          (obj["constrained_dof"][6] <> 0) then
        P := P + subs(obj["generic_support_reactions"], sol,
            TrussMe:-ComputeSpringEnergy(-MZ, (x -> obj["stiffness"](x)[6]))
          );
      end if;
    elif TrussMe:-IsCompliantJoint(obj) then
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
        P := P + subs(sol,
            TrussMe:-ComputeSpringEnergy(FJX, (x -> obj["stiffness"](x)[1]))
          );
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
        P := P + subs(sol,
            TrussMe:-ComputeSpringEnergy(FJY, (x -> obj["stiffness"](x)[2]))
          );
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
        P := P + subs(sol,
            TrussMe:-ComputeSpringEnergy(FJZ, (x -> obj["stiffness"](x)[3]))
          );
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
        P := P + subs(sol,
            TrussMe:-ComputeSpringEnergy(MJX, (x -> obj["stiffness"](x)[4]))
          );
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
        P := P + subs(sol,
            TrussMe:-ComputeSpringEnergy(MJY, (x -> obj["stiffness"](x)[5]))
          );
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
        P := P + subs(sol,
            TrussMe:-ComputeSpringEnergy(MJZ, (x -> obj["stiffness"](x)[6]))
          );
      end if;
    end if;
  end do;

  if _nresults = 1 then
    return P;
  else
    return P, uveils;
  end if;
end proc: # ComputePotentialEnergy

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
  }, $)

  description "Solve the isostatic structure equilibrium equation system given "
    "the structure objects <objs>, the external actions <exts> and the variables "
    "<vars> to solve.";

  local iso_eq, iso_eq_tmp, exts_comps, x, obj, ext, active_ext, iso_sol, A, B,
    iso_vars, veils;

  # Compute structure equations
  if (m_VerboseMode > 0) then
    printf("TrussMe:-IsostaticSolver(...): computing the equilibrium equations...\n");
  end if;
  iso_eq := [];
  for obj in objs do
    active_ext := {};
    for ext in exts do
      if (ext["target"] = obj["name"]) then
        active_ext := active_ext union {eval(ext)};
      end if;
    end do;
    iso_eq := iso_eq union TrussMe:-NewtonEuler(active_ext, obj);
    # Add joints and supports constraint equations
    if TrussMe:-IsSupport(obj) or TrussMe:-IsJoint(obj) then
      iso_eq := iso_eq union obj["constraint_loads"];
    end if;
  end do;

  # Remove NULL equations
  iso_eq := remove(x -> x = 0, TrussMe:-Simplify(iso_eq));

  # Remove equation related to rigid body motions
  # FIXME: for qloads should be different (this check can be skipped)
  iso_eq_tmp := iso_eq;
  exts_comps := map(
    x -> op(x["components"]), TrussMe:-GetObjsByType({FORCE, MOMENT}, exts)
  );
  #iso_eq := remove(x -> (
  #  member(0., (abs~(exts_comps) -~ abs(x)) *~ 1.) and (not member(0., (abs~(vars) -~ abs(x)) *~ 1.))
  #), iso_eq_tmp);
  iso_eq := remove(x -> (
    member(0., (abs~(exts_comps) -~ abs(x)) *~ 1.) and (not has(vars, indets(x)))
  ), iso_eq_tmp);
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

  if (m_VerboseMode > 0) then
    printf("TrussMe:-IsostaticSolver(...): computing the equilibrium equations... DONE\n");
    if (m_VerboseMode > 1) then
      printf("TrussMe:-IsostaticSolver(...): structure equilibrium equations:\n");
      print(<op(iso_eq)>);
      printf("TrussMe:-IsostaticSolver(...): structure unknown variables:\n");
      print(iso_vars);
    end if;
  end if;

  # Check for implicit solution flag
  if implicit then
    iso_sol := [];
  else
    if nops(iso_eq) = nops(iso_vars) then
      if (m_VerboseMode > 0) then
        if (m_VerboseMode > 1) then
          # Matrix form
          A, B := LinearAlgebra:-GenerateMatrix(iso_eq, iso_vars);
          A := Matrix(A, storage = sparse);
          printf("TrussMe:-IsostaticSolver(...): matrix visualization of the "
            "linear system:\n");
          print(plots:-sparsematrixplot(A, matrixview));
        end if;
        printf("TrussMe:-IsostaticSolver(...): computing the structure reaction "
          "forces...\n");
      end if;
      # Solve structure equations (LinearSolver)
      if m_LinearSolver = "LAST" then
        iso_sol, veils := TrussMe:-LinearSolver(iso_eq, iso_vars);
      elif m_LinearSolver = "Maple" then
        iso_sol := op(solve(iso_eq, iso_vars)) assuming real;
        veils := [];
      else
        error "TrussMe:-IsostaticSolver(...): invalid linear solver.";
      end if;

    else

      if m_WarningMode then
        WARNING("TrussMe:-IsostaticSolver(...): the system of equations is not "
          "consistent, trying to solve the system of equations anyway without "
          "LinearSolver");
      end if;
      # Solve structure equations (solve)
      iso_sol := op(solve(iso_eq, iso_vars)) assuming real;
      veils := [];
    end if;

    if (iso_sol = NULL) then
      error "TrussMe:-IsostaticSolver(...): isostatic solution not found.";
    end if;

    if (m_VerboseMode > 0) then
      printf("TrussMe:-IsostaticSolver(...): computing the structure reaction "
        "forces... DONE\n");
    end if;
  end if;

  if (_nresults = 4) then
    return iso_sol, veils, iso_eq, iso_vars;
  else
    return iso_sol, veils;
  end
end proc: # IsostaticSolver

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputeInternalActions := proc(
  objs::{
    list({BEAM, ROD}),
    set( {BEAM, ROD})
  },
  exts::{
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set( {FORCE, MOMENT, QFORCE, QMOMENT})
  },
  sol::{list, set},
  $)

  description "Programmatic computation of internal actions for structure "
    "objects <objs> with given external actions <exts> and structure solution "
    "<sol>.";

  local i, j, active_ext, subs_ext;

  # Substitute structure solution into loads
  subs_ext := map(convert, map2(subs, sol, map(op, exts)), table);

  for i from 1 to nops(objs) do
    # Extract active loads
    active_ext := {};
    for j from 1 to nops(subs_ext) do
      if (subs_ext[j]["target"] = objs[i]["name"]) then
        active_ext := active_ext union {subs_ext[j]};
      end if;
    end do;
    # Compute internal actions
    TrussMe:-InternalActions(objs[i], active_ext);
  end do;
end proc: # ComputeInternalActions

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export InternalActions := proc(
  obj::{BEAM, ROD},
  exts::{
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set({FORCE, MOMENT, QFORCE, QMOMENT})
  }, $)

  description "Programmatic computation of internal actions for structure "
    "objects <objs> with given external actions <exts> and structure solution "
    "<sol>. The function return the internal actions as function of the axial "
    "coordinate 'x'.";

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
    "TrussMe:-InternalActions(...): updating %s %s's internal actions...\n",
    obj["type"], obj["name"]
  );
  end if;

  obj["internal_actions"] := ia;

  if (m_VerboseMode > 1) then
    printf(
      "TrussMe:-InternalActions(...): updating %s %s's internal actions... DONE\n",
      obj["type"], obj["name"]
    );
  end if;

  return NULL;
end proc: # InternalActions

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputeSpringDisplacement := proc(
  spring_load::algebraic,
  spring_stiffness::procedure,
  $)::anything; # FIXME: anything is not correct

  description "Compute the displacement of a spring given the load <spring_load> "
  "and spring stiffness <stiffness>.";

  local x, out;

  # Physics:-Assume(spring_load * Dx >= 0);

  # This works even for negative Dx
  out := solve(
    spring_load = TrussMe:-Simplify(integrate(spring_stiffness(x), x = 0..Dx)),
    Dx,
    useassumptions = true,
    dropmultiplicity = true
    ) assuming real;

  # FIXME: deal with multiple solutions

  ##print("DISP EQ: ", spring_load = TrussMe:-Simplify(integrate(spring_stiffness(x), x = 0..Dx)));

  #SolveTools:-Engine({
  #  spring_load = TrussMe:-Simplify(integrate(spring_stiffness(x), x = 0..Dx))}, {Dx}, explicit):

  ##print("DISP SOL: ", %);
  #out := subs(TrussMe:-Simplify(remove(x-> has(x,I), op~(%))), Dx);
  # FIXME: deal with multiple solutions
  ##print("DISP SOL OUT: ", %);

  # Remove weird list notation inside a picewise function
  # FIXME: This works only if out is a single piece piecewise function
  if type(out, piecewise) then
    out := piecewise(op(map(
      x -> `if`(type(x, list), op(x), x), convert(out, list)
    )));
  end if;
  return out;
end proc: # ComputeSpringDisplacement

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputeSpringEnergy := proc(
  spring_load::algebraic,
  spring_stiffness::procedure,
  $)::algebraic;

  description "Compute the potential energy of a spring given the load "
  "<spring_load> and spring stiffness <stiffness>.";

  local disp;

  disp := TrussMe:-ComputeSpringDisplacement(spring_load, spring_stiffness);
  return TrussMe:-Simplify(
    integrate(integrate(spring_stiffness(x), x), x = 0..disp)
  );
end proc: # ComputeSpringEnergy

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputeSupportDisplacements := proc(
  obj::SUPPORT,
  $)

  description "Compute the displacements of the support <obj> from its "
    "support reactions.";

  local sup_disp, disp_vec, i, disp, x, sup_reac;

  sup_disp := [];
  disp_vec := [tx, ty, tz, rx, ry, rz];
  sup_reac := [FX, FY, FZ, MX, MY, MZ];

  for i from 1 to 6 do
    if (obj["constrained_dof"][i] = 1) and
        member(sup_reac[i], map(lhs, obj["generic_support_reactions"])) then
      disp := TrussMe:-ComputeSpringDisplacement(
        subs(obj["support_reactions"], - sup_reac[i]),
        (x -> obj["stiffness"](x)[i])
      );
      sup_disp := sup_disp union [disp_vec[i] = disp];
    end if;
  end do;

  if (m_VerboseMode > 1) then
    printf(
      "TrussMe:-ComputeSupportDisplacements(...): updating %s %s's displacements...\n",
      obj["type"], obj["name"]
    );
  end if;

  obj["displacements"] := sup_disp;

  if (m_VerboseMode > 1) then
      "TrussMe:-ComputeSupportDisplacements(...): updating %s %s's displacements... "
      "DONE\n", obj["type"], obj["name"]
  end if;
  return NULL;
end proc: # ComputeSupportDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputeJointDisplacements := proc(
  obj::JOINT,
  sol::{list, set},
  $)

  description "Compute the displacements of the joint <obj> given the solution "
    "<sol>.";

  local jnt_disp, disp_vec, i, disp, x, jnt_load, f;

  jnt_disp := [];
  disp_vec := ['tx', 'ty', 'tz', 'rx', 'ry', 'rz'];
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
      disp := TrussMe:-ComputeSpringDisplacement(subs(sol, jnt_load[i]),
        (x -> obj["stiffness"](x)[i]));
      jnt_disp := jnt_disp union [disp_vec[i] = disp];
    end if;
  end do;

  if (m_VerboseMode > 1) then
    printf(
      "TrussMe:-ComputeJointDisplacements(...): updating %s %s's displacements...\n",
      obj["type"], obj["name"]
    );
  end if;

  obj["displacements"] := jnt_disp;

  if (m_VerboseMode > 1) then
    printf(
      "TrussMe:-ComputeJointDisplacements(...): updating %s %s's displacements... "
      "DONE\n", obj["type"], obj["name"]
    );
  end if;

return NULL;
end proc: # ComputeJointDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputeDisplacements := proc(
  objs::{
    list({BEAM, ROD, SUPPORT, JOINT}),
    set({BEAM, ROD, SUPPORT, JOINT})
  },
  exts::{
    list({FORCE, MOMENT, QFORCE, QMOMENT}),
    set({FORCE, MOMENT, QFORCE, QMOMENT})
  },
  sol::{list, set},
  {
    timoshenko_beam::boolean := false
  }, $)

  description "Compute the structure displacements and rotations given the "
    "structure objects <objs>, the exteranl forces <exts>, the solution <sol>, "
    "and the Timoshenko beam flag <timoshenko_beam>.";

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
    elif TrussMe:-IsCompliantSupport(obj) then
      # Compute displacements
      TrussMe:-ComputeSupportDisplacements(obj);

    # Joint
    elif TrussMe:-IsCompliantJoint(obj) then
      # Compute displacements
      TrussMe:-ComputeJointDisplacements(obj, sol);
    end if;
  end do;
end proc: # ComputeDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputePunctualDisplacement := proc(
  struct::STRUCTURE,
  objs::{
    list({BEAM, ROD, RIGID_BODY, SUPPORT, JOINT})
  },
  coords::list,
  directions::list,
  RFs::list,
  {
    timoshenko_beam::boolean := false,
    unveil_results::boolean  := true
  }, $)

  description "Compute the structure <struct> punctual displacements of the "
    "object <obj> at the coordinates <coords> in the directions <directions>. "
    "The directions are defined in the reference frame <RFs>. Optional argument "
    "are: <timoshenko_beam> boolean flag to use Timoshenko beam model, and "
    "<unveil_results> boolean flag to unveil the results.";

  local out, struct_copy, obj, objs_names, dummy_loads, subs_obj, obj_coords,
    obj_targets, x, subs_null_dummy, disp, i, j, d_coords, sw_tmp, out_veils;

  # Substitute -1 entries of coords with the corresponding object length
  d_coords := [seq(`if`(coords[i] = -1, objs[i]["length"], coords[i]), i = 1..nops(coords))];

  # Set module local variable m_KeepVeiled
  m_KeepVeiled := not unveil_results;

  # Create a copy of the structure
  struct_copy := TrussMe:-CopyStructure(struct);
  TrussMe:-CleanStructure(struct_copy);
  # Get objects names
  objs_names := TrussMe:-GetNames(objs);

  # Disable warnings temporarily
  sw_tmp := m_WarningMode;
  m_WarningMode := false;
  # Replace rods with beams to be able to compute the displacements in all directions
  for obj in map(TrussMe:-GetObjByName, objs_names, struct_copy["objects"]) do
    if TrussMe:-IsRod(obj) then
      subs_obj := TrussMe:-MakeBeam(
        obj["name"],
        obj["length"],
        obj["frame"],
        parse("area") = obj["area"],
        parse("material") = obj["material"]
      );
      # Remove load on unconstrained direction
      subs_obj["admissible_loads"] := [1, 1, 1, 0, 1, 1];
      # Replace object in struct_copy
      struct_copy["objects"] := remove(
        x -> x["name"] = obj["name"], struct_copy["objects"]
      );
      struct_copy["objects"] := struct_copy["objects"] union {eval(subs_obj)};
    end if;
  end do;
  # Update struct_copy supports and joints for the new objects
  for obj in struct_copy["objects"] do
    if TrussMe:-IsSupport(obj) then
      # Re-make support to generate new loads and constraint compliant with
      # substituted objects (joint is made because earth is already in the list
      # of targets)
      obj_targets := map(TrussMe:-GetObjByName, remove(
        x -> x = m_earth["name"], obj["targets"]
      ), struct_copy["objects"]);
      obj_coords  := obj["coordinates"][2..-1];
      subs_obj := TrussMe:-MakeSupport(
        obj["name"],
        obj["constrained_dof"],
        obj_targets,
        obj_coords,
        obj["frame"],
        parse("stiffness") = obj["stiffness"]
      );
      # Replace object in struct_copy
      struct_copy["objects"] := remove(
        x -> x["name"] = obj["name"], struct_copy["objects"]
      );
      struct_copy["objects"] := struct_copy["objects"] union {eval(subs_obj)};
    elif TrussMe:-IsJoint(obj) then
      # Re-Make joint to generate new loads and constraint compliant with substituted objects
      obj_targets := map(TrussMe:-GetObjByName, obj["targets"], struct_copy["objects"]);
      subs_obj := TrussMe:-MakeJoint(
        obj["name"],
        obj["constrained_dof"],
        obj_targets,
        obj["coordinates"],
        obj["frame"],
        parse("stiffness") = obj["stiffness"]
      );
      # Replace object in struct_copy
      struct_copy["objects"] := remove(
        x -> x["name"] = obj["name"], struct_copy["objects"]
      );
      struct_copy["objects"] := struct_copy["objects"] union {eval(subs_obj)};
    end if;
  end do;

  # Create dummy loads in the directions of interest
  subs_null_dummy := [];
  for i from 1 to nops(objs_names) do
    dummy_loads := eval~([
      `if`(directions[i, 1] = 1, TrussMe:-MakeForce( [dFx_||i,0,0], d_coords[i], TrussMe:-GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL),
      `if`(directions[i, 2] = 1, TrussMe:-MakeForce( [0,dFy_||i,0], d_coords[i], TrussMe:-GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL),
      `if`(directions[i, 3] = 1, TrussMe:-MakeForce( [0,0,dFz_||i], d_coords[i], TrussMe:-GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL),
      `if`(directions[i, 4] = 1, TrussMe:-MakeMoment([dMx_||i,0,0], d_coords[i], TrussMe:-GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL),
      `if`(directions[i, 5] = 1, TrussMe:-MakeMoment([0,dMy_||i,0], d_coords[i], TrussMe:-GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL),
      `if`(directions[i, 6] = 1, TrussMe:-MakeMoment([0,0,dMz_||i], d_coords[i], TrussMe:-GetObjByName(objs_names[i], struct_copy["objects"]), RFs[i]), NULL)
    ]);

    # Null dummy loads substitution list
    subs_null_dummy := subs_null_dummy union (
      [dFx_||i, dFy_||i, dFz_||i, dMx_||i, dMy_||i, dMz_||i] =~ [0, 0, 0, 0, 0, 0]
    );

    # Add dummy loads to the structure copy
    struct_copy["external_actions"] := struct_copy["external_actions"] union dummy_loads;
  end do;

  # Solve the structure copy
  TrussMe:-SolveStructure(
    struct_copy,
    parse("compute_internal_actions") = false,
    parse("compute_displacements")    = false,
    parse("compute_potential_energy") = true,
    parse("timoshenko_beam")          = timoshenko_beam,
    parse("implicit")                 = false,
    parse("unveil_results")           = unveil_results,
    parse("dummy_vars")               = lhs~(subs_null_dummy)
    );

  # Reset warning mode
  m_WarningMode := sw_tmp;

  # Compute punctual displacements
  if (m_VerboseMode > 0) then
    printf("TrussMe:-ComputePunctualDisplacement: Computing punctual "
      "displacements... \n");
  end if;

  out := [];
  for i from 1 to nops(objs_names) do
    disp := [];
    if directions[i,1] = 1 then
      disp := disp union [ux = subs(subs_null_dummy, TrussMe:-Diff(
          subs(subs_null_dummy[[2,3,4,5,6]], struct_copy["potential_energy"]),
          dFx_||i,
          parse("veils") = subs(subs_null_dummy[[2,3,4,5,6]], struct_copy["veils"])
        ))];
    end if;
    if directions[i,2] = 1 then
      disp := disp union [uy = subs(subs_null_dummy, TrussMe:-Diff(
          subs(subs_null_dummy[[1,3,4,5,6]], struct_copy["potential_energy"]),
          dFy_||i,
          parse("veils") = subs(subs_null_dummy[[1,3,4,5,6]], struct_copy["veils"])
        ))];
    end if;
    if directions[i,3] = 1 then
      disp := disp union [uz = subs(subs_null_dummy, TrussMe:-Diff(
          subs(subs_null_dummy[[1,2,4,5,6]], struct_copy["potential_energy"]),
          dFz_||i,
          parse("veils") = subs(subs_null_dummy[[1,2,4,5,6]], struct_copy["veils"])
        ))];
    end if;
    if directions[i,4] = 1 then
      disp := disp union [rx = subs(subs_null_dummy, TrussMe:-Diff(
          subs(subs_null_dummy[[1,2,3,5,6]], struct_copy["potential_energy"]),
          dMx_||i,
          parse("veils") = subs(subs_null_dummy[[1,2,3,5,6]], struct_copy["veils"])
        ))];
    end if;
    if directions[i,5] = 1 then
      disp := disp union [ry = subs(subs_null_dummy, TrussMe:-Diff(
          subs(subs_null_dummy[[1,2,3,4,6]], struct_copy["potential_energy"]),
          dMy_||i,
          parse("veils") = subs(subs_null_dummy[[1,2,3,4,6]], struct_copy["veils"])
        ))];
    end if;
    if directions[i,6] = 1 then
      disp := disp union [rz = subs(subs_null_dummy, TrussMe:-Diff(
          subs(subs_null_dummy[[1,2,3,4,5]], struct_copy["potential_energy"]),
          dMz_||i,
          parse("veils") = subs(subs_null_dummy[[1,2,3,4,5]], struct_copy["veils"])
        ))];
    end if;
    # Select displacement relative to desired directions and add to the list of displacements
    out := out union [disp];
  end do;

  if (m_VerboseMode > 0) then
    printf("TrussMe:-ComputePunctualDisplacement: Computing punctual "
      "displacements... DONE.\n");
  end if;

  # Compose output veils
  # FIXME: only veils which are present in 'out' can be returned
  # (with their recursive dependencies), They should not contain dummy variables
  out_veils := subs(subs_null_dummy, struct_copy["veils"]);

  # Simplify output
  out := TrussMe:-Simplify(out);

  if _nresults = 1 then
    return out;
  else
    return out, out_veils;
  end if;
end proc: # ComputePunctualDisplacement

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export ComputeObjectFrameDisplacements := proc(
  struct::STRUCTURE,
  {
    timoshenko_beam::boolean := false,
    unveil_results::boolean  := true
  }, $)::list;

description "Compute the total displacements of the structure <struct> with "
  "optional <timoshenko_beam> and <unveil_results> flags.";

  local i, RF_nt, nx, ny, nz, theta, subs_n, subs_t, disp, veils, objs;

  # Compute punctual displacements at origin for all the structure objects
  # FIXME: this is overkill, we should compute only the displacements in the
  # direction of joints dof (and not even all of them)
  disp, veils := TrussMe:-ComputePunctualDisplacement(
    struct,
    convert(struct["objects"],list),
    [seq([0, 0, 0], i = 1..nops(struct["objects"]))],
    [seq([1, 1, 1, 1, 1, 1], i = 1..nops(struct["objects"]))],
    map(x-> x["frame"], convert(struct["objects"], list)),
    parse("timoshenko_beam") = timoshenko_beam,
    parse("unveil_results")  = unveil_results
    ):

  # Update objects frame displacements
  for i from 1 to nops(struct["objects"]) do
    # FIXME: not all displacement are necessarily computed rotations about a
    # generic axis
    RF_nt  := Matrix(4, 4, [[
        -nx^2*cos(theta) + nx^2 + cos(theta),
        -nx*ny*cos(theta) - sin(theta)*nz + nx*ny,
        -nx*nz*cos(theta) + sin(theta)*ny + nx*nz, 0
      ], [
        -nx*ny*cos(theta) + sin(theta)*nz + nx*ny,
        (-nx^2 + nz^2 + 1)*cos(theta) + nx^2 - nz^2,
        -sin(theta)*nx - ny*nz*(cos(theta) - 1), 0
      ], [
        -nx*nz*cos(theta) - sin(theta)*ny + nx*nz,
        sin(theta)*nx - ny*nz*(cos(theta) - 1),
        -cos(theta)*(-2*nx^2 + nz^2) - 2*nx^2 + nz^2 + 1, 0
      ], [
        0, 0, 0, 1
    ]]);
    if subs(disp[i], TrussMe:-Norm2([rx, ry, rz])) <> 0. and
       subs(disp[i], TrussMe:-Norm2([rx, ry, rz])) <> 0 then
      subs_n := [nx, ny, nz] =~ subs(
        disp[i], [rx, ry, rz] /~ TrussMe:-Norm2([rx, ry, rz])
      );
    else
      subs_n := [nx, ny, nz] =~ [0, 0, 1];
    end if;
    subs_t := theta = subs(disp[i], TrussMe:-Norm2([rx, ry, rz]));
    struct["objects"][i]["frame_displacements"] :=
      TrussMe:-Translate(
        op(subs(disp[i], [ux, uy, uz]))).subs(subs_n, subs_t, RF_nt
      );
  end do;

  return veils;
end proc: # ComputeObjectFrameDisplacements

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

export LinearSolver := proc(
  eqns::{list, set},
  vars::{list, set},
  $)

  description "Solve the linear system of equations <eqns> for the variables "
    "<vars>.";

  local T, sol, sol_tmp, veils, A, B;

  print(m_StoredData);

  # Set veiling_label increasing index
  m_VeilingIdx := m_VeilingIdx + 1;
  m_LEM:-SetVeilingLabel(m_LEM, cat(m_VeilingLabel, m_VeilingIdx));

  # Matrix form of the linear system
  A, B := LinearAlgebra:-GenerateMatrix(eqns, vars);

  if (has(A, vars)) or (has(B, vars)) then
    error "TrussMe:-LinearSolver(...): the system is not linear in the "
      "variables %1.", vars;
  end if;

  if (m_VerboseMode > 0) then
    printf("TrussMe:-LinearSolver(...): performing LU decomposition...\n");
  end if;

  # LU decomposition
  if (m_StoredData <> NULL) then
    m_LAST:-SetStoredData(m_LAST, m_StoredData);
  end if;
  m_LAST:-LU(m_LAST, A);

  if (m_VerboseMode > 0) then
    printf("TrussMe:-LinearSolver(...): performing LU decomposition... DONE\n");
    printf("TrussMe:-LinearSolver(...): solving linear system...\n");
  end if;

  # Solve linear system
  sol_tmp := m_LAST:-SolveLinearSystem(m_LAST, B);

  if (m_VerboseMode > 0) then
    printf("TrussMe:-LinearSolver(...): solving linear system... DONE\n");
    if (m_VerboseMode > 1) then
      print("sol_tmp", sol_tmp);
    end if;
  end if;

  # Substitute veils to solution
  if m_KeepVeiled then
    # Remove indexed type from veils
    m_LEM:-VeilList(m_LEM);
    lhs~(%) =~ map2(op, 0, lhs~(%)) ||~ __ ||~ (op~(lhs~(%)));
    # Substitutution
    sol := convert(vars =~ subs(%, sol_tmp), list);
    veils := subs(%%, %%%);
    # Update StoredData
    subs(op(m_StoredData), veils);
    m_StoredData := m_StoredData union (lhs~(%) =~ subs(op(ListTools:-Reverse(%)), rhs~(%)));
  else
    if (m_VerboseMode > 0) then
      printf("TrussMe:-LinearSolver(...): substituting veils...\n");
    end if;
    # Substitutution
    sol := convert(vars =~ m_LEM:-UnVeil(m_LEM, sol_tmp), list);
    veils := [];
    if (m_VerboseMode > 0) then
      printf("TrussMe:-LinearSolver(...): substituting veils... DONE\n");
    end if;
  end if;
  m_LEM:-ForgetVeil(m_LEM);

  if _nresults = 1 then
    return sol;
  elif _nresults = 2 then
    return sol, veils;
  end if;
  return TrussMe:-Simplify(sol), veils;
end proc: # LinearSolver

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
