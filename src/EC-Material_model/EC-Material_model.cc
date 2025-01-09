

#include "EC-Material_model.h"

#include <aspect/utilities.h>
#include <aspect/newton.h>

#include <deal.II/numerics/fe_field_function.h>


using namespace std;
using namespace aspect::MaterialModel;

template <int dim> void EC_Material_model<dim>::evaluate ( const MaterialModel::MaterialModelInputs<dim> & in, MaterialModel::MaterialModelOutputs<dim> & out ) const
{
  // Inicializace a pripravy
  AssertThrow ( this->introspection().compositional_name_exists("mantle_depletion") == true,
    ExcMessage ( "mantle_depletion composition field not detected!" ) );
  EquationOfStateOutputs<dim> eos_outputs (1);
  const int mantle_depletion_indx = this->introspection().compositional_index_for_name("mantle_depletion");

  // Pripravit dodatecne outputy a reactionRateOutputy
  create_additional_named_outputs ( out );
  EC_Additional_outputs<dim> * AO = out.template get_additional_output<EC_Additional_outputs<dim>>();
  AssertThrow ( AO != nullptr, ExcMessage("") );
  ReactionRateOutputs<dim> * reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim>>();
  AssertThrow ( this->get_parameters().use_operator_splitting == false || reaction_rate_out != nullptr, 
    ExcMessage ( string("this->get_parameters().use_operator_splitting = ") + (this->get_parameters().use_operator_splitting == true ? "true" : "false") + ",\nreaction_rate_out = " + (reaction_rate_out == nullptr ? "nullptr" : "not nullptr") ) );

  // Nalezeni stareho F
  vector<double> Fold ( in.n_evaluation_points(), 0.0 );
  if ( par_melt.calculate_depletion == true && in.current_cell.state() == IteratorState::valid && this->get_timestep_number() > 0 )
  {
    Functions::FEFieldFunction<dim, LinearAlgebra::BlockVector> fe_value ( this->get_dof_handler(), this->get_old_solution(), this->get_mapping() );
    fe_value.set_active_cell ( in.current_cell );
    fe_value.value_list ( in.position, Fold, this->introspection().component_indices.compositional_fields[mantle_depletion_indx] );
  }

  

  for ( unsigned int i = 0; i < in.n_evaluation_points(); i++ )
  {
    equation_of_state         .evaluate ( in, i, eos_outputs );
    const double pressure     = max ( in.pressure[i], 0.0 );
    const double temperature  = in.temperature[i];

    // Constanty z clanku Katz 2003
    const double cA [3] = { 1085.7, 132.9e-9, -5.1e-18 };
    const double cC [3] = { 1780.0, 45.0e-9,  -2.0e-18 };

    // Solidus a Liquidus z clanku Katz 2003, prevedeno z Celsia na Kelviny
    const double CtoK = 273.14;
    const double Sol  = cA[0] + cA[1]*pressure + cA[2]*pow ( pressure, 2.0 ) + CtoK;
    const double Liq  = cC[0] + cC[1]*pressure + cC[2]*pow ( pressure, 2.0 ) + CtoK;

    // Napocteni parametricke teploty dle McKenzie 1988
    const double depth          = this->get_geometry_model().depth(in.position[i]);
    const double alpha          = eos_outputs.compressibilities[0];
    const double specific_heat  = eos_outputs.specific_heat_capacities[0];
    const double g              = this->get_gravity_model().gravity_vector(in.position[i]).norm();
    const double Ta             = temperature + depth*alpha*g*temperature/specific_heat;
    const double Tss            = min ( 0.5, max ( ( Ta - (Sol + Liq)/2.0 ) / (Liq - Sol), -0.5 ) );

    // Napocteni hodnoty meltu F dle McKenzie 1988
    const double a = 0.5;
    const double b = 0.25;
    const double c = 0.4259;
    const double d = 2.988;
    const double F = a + Tss + ( Tss*Tss - b ) * ( c + d*Tss );

    // Nalezeni soucasneho F
    const double Fnow = max ( F, Fold[i] );

    // Napocteni teploty pro viskozitu dle Capitanio 2024 (odecteno fazove teplo?)
    const double Tstar = Ta - Fnow*260.0;

    AssertThrow ( AO->liquidus.size() == in.n_evaluation_points(),
      ExcMessage ( "AO->liquidus.size() = " + to_string(AO->liquidus.size()) + ";\nin.n_evaluation_points() = " + to_string(in.n_evaluation_points()) + ";\nout.n_evaluation_points() = " + to_string(out.n_evaluation_points()) ) );
    AO->liquidus[i] = Liq;
    AO->solidus[i] = Sol;
    AO->Ta[i] = Ta;
    AO->Tss[i] = Tss;
    AO->Finter[i] = F;
    AO->Tstar[i] = Tstar;

    // Napocitani viskozity dle modelu
    double eta_eff = par_viscosity.eta_0;
    if ( par_viscosity.use_diffusion == true )
      eta_eff *= par_viscosity.diffusion.A * exp ( par_viscosity.diffusion.E_a/R/Tstar );
    if ( par_viscosity.use_melt_depletion == true )
      eta_eff *= pow ( 1.0-Fnow, ( par_viscosity.depletion.k - 1.0 ) / par_viscosity.depletion.k );
    if ( par_viscosity.use_druckerprager == true )
    {
      const double edot_ii = max(sqrt(fabs(second_invariant(deviator(in.strain_rate[i])))), par_viscosity.druckerprager.em );
      const double yield =  ( par_viscosity.druckerprager.use_melt_depletion ) == true ?
                            ( par_viscosity.druckerprager.C + Fnow*par_viscosity.druckerprager.delta_C_deplt + ( par_viscosity.druckerprager.mu + Fnow*par_viscosity.druckerprager.delta_mu_deplt ) * pressure ) :
                            ( par_viscosity.druckerprager.C + par_viscosity.druckerprager.mu * pressure );
      eta_eff = min ( eta_eff, yield/2.0/edot_ii );
    }
    out.viscosities[i] = max ( min ( eta_eff, par_viscosity.limits.eta_max ), par_viscosity.limits.eta_min );

    // Napocitani reaction rate outputu
    if ( par_melt.calculate_depletion == true )
    {
      const double deltaF = ( ( par_melt.use_depth_regeneration == true && depth > par_melt.regeneration_depth ) ? ( -Fold[i] ) : max ( 0.0, F - Fold[i] ) );
      if ( this->get_parameters().use_operator_splitting == false )
        out.reaction_terms[i][mantle_depletion_indx] = deltaF;
      else
      {
        reaction_rate_out->reaction_rates[i][mantle_depletion_indx] = deltaF / this->get_timestep();
        out.reaction_terms[i][mantle_depletion_indx] = 0.0;
      }
    }
    else
      out.reaction_terms[i][mantle_depletion_indx] = 0.0;
    
    // Napocitani hustoty dle modelu
    double rho_eff = eos_outputs.densities[0];
    if ( par_density.use_melt_depletion == true )
      rho_eff += Fnow * par_density.depletion.delta_rho;
    out.densities[i] = rho_eff;

    // Ostatni
    out.thermal_expansion_coefficients[i] = eos_outputs.thermal_expansion_coefficients[0];
    out.specific_heat[i] = specific_heat;
    out.thermal_conductivities[i] = thermal_conductivity;
    out.compressibilities[i] = alpha;
    out.entropy_derivative_pressure[i] = eos_outputs.entropy_derivative_pressure[0];
    out.entropy_derivative_temperature[i] = eos_outputs.entropy_derivative_temperature[0];
  }
}

template <int dim> bool EC_Material_model<dim>::is_compressible () const
{
  return equation_of_state.is_compressible();
}

template <int dim> void EC_Material_model<dim>::declare_parameters ( ParameterHandler & prm )
{
  prm.enter_subsection ( "Material model" );
  {
    prm.enter_subsection ( "Enhanced capitanio" );
    {
      EquationOfState::LinearizedIncompressible<dim>::declare_parameters (prm);
      prm.declare_entry ( "Thermal conductivity", "3.96", Patterns::Double(0.0),
                          "The value of the thermal conductivity $k = \\kappa c_p \\alpha$, where $\\kappa$ is thermal diffusivity. Units: \\si{\\watt\\per\\meter\\per\\kelvin}." );
      prm.enter_subsection ( "Melt options" );
      {
        prm.declare_entry ( "Calculate depletion", "true", Patterns::Bool(),
                            "Whether amount of depletion should be calculated and advected." );
        prm.declare_entry ( "Use depth regeneration", "false", Patterns::Bool(),
                            "Whether the depleted mantle material should regenerate when beneath regeneration depth. If set to true, amount of depletion is set to zero for all mantle material in depths greater than regeneration depth." );
        prm.declare_entry ( "Regeneration depth", "300.0e3", Patterns::Double (),
                            "Regeneration depth, check Use depth regeneration parameter. Units: \\si{\\meter}." );
      }
      prm.leave_subsection ();
      prm.enter_subsection ( "Viscosity" );
      {
        const string available_viscosity_models = "none|drucker-prager|melt depletion|diffusion";
        prm.declare_entry ( "List of viscosity models", "drucker-prager, melt depletion, diffusion", Patterns::MultipleSelection ( available_viscosity_models ),
                            "The list of models used when calculating viscosity. Models available: " + available_viscosity_models );
        prm.declare_entry ( "Reference viscosity", "1.0e20", Patterns::Double (0.0),
                            "The value of the reference viscosity $\\eta_0$. Units: \\si{\\pascal\\second}." );
        prm.enter_subsection ( "Diffusion creep" );
        {
          prm.declare_entry ( "Prefactor", "9.1963e-9", Patterns::Double (0.0),
                              "The value of the prefactor $A$." );
          prm.declare_entry ( "Activation energy", "2.4e5", Patterns::Double (0.0),
                              "The value of the activation energy $E_a$. Units: \\si{\\joules\\per\\mol}." );
        }
        prm.leave_subsection ();
        prm.enter_subsection ( "Viscosity limits" );
        {
          prm.declare_entry ( "Minimum viscosity", "1.0e20", Patterns::Double (0.0),
                              "The value of minimum viscosity $\\eta_{min}$. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ( "Maximum viscosity", "1.0e25", Patterns::Double (0.0),
                              "The value of maximum viscosity $\\eta_{max}$. Units: \\si{\\pascal\\second}.");
        }
        prm.leave_subsection ();
        prm.enter_subsection ( "Melt depletion" );
        {
          prm.declare_entry ( "k water", "2.1e-2", Patterns::Double (0.0),
                              "The value of k. Units: Non-dimensional.");
        }
        prm.leave_subsection ();
        prm.enter_subsection ( "Drucker-prager" );
        {
          prm.declare_entry ( "Cohesion", "10.0e6", Patterns::Double (0.0),
                            "The value of cohesion. Units: \\si{\\mega\\pascal}");
          prm.declare_entry ( "Friction", "0.013", Patterns::Double (0.0),
                            "The value of friction coefficient. Units: No units");
          prm.declare_entry ( "Minimum strain rate", "1.0e-16", Patterns::Double (),
                            "The value of minimum strain rate. Units: \\si{\\second^{-1}}");
          prm.declare_entry ( "Use melt depletion", "false", Patterns::Bool (),
                            "Whether melt depletion changes cohesion and friction." );
          prm.declare_entry ( "Depleted cohesion change", "0.0", Patterns::Double (0.0),
                            "The value of cohesion change of depleted material. Units: \\si{\\mega\\pascal}" );
          prm.declare_entry ( "Depleted friction change", "0.0", Patterns::Double (0.0),
                            "The value of friction coefficient change of depleted material. Units: No units");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
      prm.enter_subsection ( "Density" );
      {
        const string available_density_models = "none|melt depletion";
        prm.declare_entry ( "List of density models", "melt depletion", Patterns::MultipleSelection ( available_density_models ),
                            "The list of models used when calculating density. Models available: " + available_density_models );
        prm.enter_subsection ( "Melt depletion" );
        {
          prm.declare_entry ( "Density change", "-72.6", Patterns::Double (),
                              "The value of density change due to depletion of mantle material. Units: \\si{\\kilogram\\per\\meter^3}");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
    prm.leave_subsection ();
  }
  prm.leave_subsection ();
}

template <int dim> void EC_Material_model<dim>::parse_parameters ( ParameterHandler & prm )
{
  prm.enter_subsection ( "Material model" );
    prm.enter_subsection ( "Enhanced capitanio" );
      equation_of_state                                   . parse_parameters ( prm );
      thermal_conductivity                                = prm.get_double ( "Thermal conductivity" );
      prm.enter_subsection ( "Melt options" );
        par_melt.calculate_depletion                      = prm.get_bool ( "Calculate depletion" );
        par_melt.use_depth_regeneration                   = prm.get_bool ( "Use depth regeneration" );
        par_melt.regeneration_depth                       = prm.get_double ( "Regeneration depth" );
      prm.leave_subsection ();
      prm.enter_subsection ( "Viscosity" );
        vector<string> vm                                 = Utilities::split_string_list ( prm.get ( "List of viscosity models" ) );
        par_viscosity.use_druckerprager                   = find ( vm.begin(), vm.end(), "drucker-prager" ) == vm.end() ? false : true;
        par_viscosity.use_melt_depletion                  = find ( vm.begin(), vm.end(), "melt depletion" ) == vm.end() ? false : true;
        par_viscosity.use_diffusion                       = find ( vm.begin(), vm.end(), "diffusion" ) == vm.end() ? false : true;
        par_viscosity.eta_0                               = prm.get_double ( "Reference viscosity" );
        prm.enter_subsection ( "Diffusion creep" );
          par_viscosity.diffusion.A                       = prm.get_double ( "Prefactor" );
          par_viscosity.diffusion.E_a                     = prm.get_double ( "Activation energy" );
        prm.leave_subsection ();
        prm.enter_subsection ( "Viscosity limits" );
          par_viscosity.limits.eta_min                    = prm.get_double ( "Minimum viscosity" );
          par_viscosity.limits.eta_max                    = prm.get_double ( "Maximum viscosity" );
        prm.leave_subsection ();
        prm.enter_subsection ( "Melt depletion" );
          par_viscosity.depletion.k                       = prm.get_double ( "k water" );
        prm.leave_subsection ();
        prm.enter_subsection ( "Drucker-prager" );
          par_viscosity.druckerprager.C                   = prm.get_double ( "Cohesion" );
          par_viscosity.druckerprager.mu                  = prm.get_double ( "Friction" );
          par_viscosity.druckerprager.em                  = prm.get_double ( "Minimum strain rate" );
          par_viscosity.druckerprager.use_melt_depletion  = prm.get_bool ( "Use melt depletion" );
          par_viscosity.druckerprager.delta_C_deplt             = prm.get_double ( "Depleted cohesion change" );
          par_viscosity.druckerprager.delta_mu_deplt            = prm.get_double ( "Depleted friction change" );
        prm.leave_subsection ();
      prm.leave_subsection ();
      prm.enter_subsection ( "Density" );
        vector<string> dm                                 = Utilities::split_string_list ( prm.get ( "List of density models" ) );
        par_density.use_melt_depletion                    = find ( dm.begin(), dm.end(), "melt depletion" ) == dm.end() ? false : true;
        prm.enter_subsection ( "Melt depletion" );
          par_density.depletion.delta_rho                 = prm.get_double ( "Density change" );
        prm.leave_subsection ();
      prm.leave_subsection ();
    prm.leave_subsection ();
  prm.leave_subsection ();

  this->model_dependence.compressibility = NonlinearDependence::none;
  this->model_dependence.specific_heat = NonlinearDependence::none;
  this->model_dependence.thermal_conductivity = NonlinearDependence::none;
  this->model_dependence.viscosity = NonlinearDependence::pressure | NonlinearDependence::temperature | NonlinearDependence::compositional_fields;
  this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::compositional_fields;
}

template <int dim> void EC_Material_model<dim>::create_additional_named_outputs ( MaterialModel::MaterialModelOutputs<dim> & out ) const
{
  if ( out.template get_additional_output<EC_Additional_outputs<dim>>() == nullptr )
  {
    const unsigned int N = out.n_evaluation_points();
    out.additional_outputs.emplace_back ( make_unique<MaterialModel::EC_Additional_outputs<dim>>(N) );
  }
  if ( this->get_parameters().use_operator_splitting == true && out.template get_additional_output<ReactionRateOutputs<dim>>() == nullptr )
  {
    const unsigned int N = out.n_evaluation_points();
    out.additional_outputs.emplace_back ( make_unique<MaterialModel::ReactionRateOutputs<dim>> ( N, this->n_compositional_fields() ) );
  }
}


namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(EC_Material_model,
                                   "enhanced capitanio",
                                   "")
  }
}
