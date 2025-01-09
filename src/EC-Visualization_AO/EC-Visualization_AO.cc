
#include "EC-Visualization_AO.h"

#include <aspect/postprocess/visualization/material_properties.h>
#include <aspect/material_model/interface.h>
#include <aspect/utilities.h>
#include <aspect/melt.h>

#include <algorithm>

using namespace aspect;
using namespace aspect::Postprocess::VisualizationPostprocessors;
using namespace std;

template <int dim> vector<string> EC_Visualization_AO<dim>::get_names () const
{
  vector<string> names;
  for ( const string & prop: property_names )
  {
    if ( prop == "solidus" )
      names.push_back("solidus");
    else if ( prop == "liquidus" )
      names.push_back("liquidus");
    else if ( prop == "T adiabatic" )
      names.push_back("T adiabatic");
    else if ( prop == "T parametric" )
      names.push_back("T parametric");
    else if ( prop == "F parametric" )
      names.push_back("F parametric");
    else if ( prop == "T star" )
      names.push_back("T star");
  }
  return names;
}

template <int dim> vector<DataComponentInterpretation::DataComponentInterpretation> EC_Visualization_AO<dim>::get_data_component_interpretation () const
{
  vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;
  for ( const string & prop: property_names )
    interpretation.push_back ( DataComponentInterpretation::component_is_scalar );
  return interpretation;
}

template <int dim> UpdateFlags EC_Visualization_AO<dim>::get_needed_update_flags () const
{
  return update_values | update_quadrature_points;
}

template <int dim> void EC_Visualization_AO<dim>::declare_parameters ( ParameterHandler & prm )
{
  // throw runtime_error ( "declare_parameters" );
  prm.enter_subsection("Postprocess");
  {
    prm.enter_subsection("Visualization");
    {
      prm.enter_subsection("Enhanced capitanio");
      {
        const string pattern_of_names = "liquidus|solidus|T adiabatic|T parametric|F parametric|T star";

        prm.declare_entry("List of properties",
                          "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma separated list of Enhanced Capitanio properties that should be "
                          "written whenever writing graphical output. "
                          "The following Enhanced Capitanio properties are available:\n\n"
                          +
                          pattern_of_names);
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();
}

template <int dim> void EC_Visualization_AO<dim>::evaluate_vector_field ( const DataPostprocessorInputs::Vector<dim> & input_data, std::vector<Vector<double>> & computed_quantities ) const
{
  const unsigned int n_quadrature_points = input_data.solution_values.size();
  Assert ( computed_quantities.size() == n_quadrature_points, ExcInternalError() );
  Assert ( input_data.solution_values[0].size() == this->introspection().n_components, ExcInternalError() );

  MaterialModel::MaterialModelInputs<dim> in ( input_data, this->introspection() );
  MaterialModel::MaterialModelOutputs<dim> out ( n_quadrature_points, this->n_compositional_fields() );

  this->get_material_model().evaluate ( in, out );
  MaterialModel::EC_Additional_outputs<dim> * AO = out.template get_additional_output<MaterialModel::EC_Additional_outputs<dim>>();

  for ( unsigned int q = 0; q < n_quadrature_points; ++q )
  {
    unsigned output_index = 0;
    for ( unsigned int i = 0; i < property_names.size(); ++i, ++output_index )
    {
      if ( property_names[i] == "solidus" )
      {
        AssertThrow ( q < AO->solidus.size(),
          ExcMessage ( "EC_Visualization_AO::evaluate_vector_field -- q >= AO->solidus.size()" ) );
        computed_quantities[q][output_index] = AO->solidus[q];
      }
      else if ( property_names[i] == "liquidus" )
      {
        AssertThrow ( q < AO->liquidus.size(),
          ExcMessage ( "EC_Visualization_AO::evaluate_vector_field -- q >= CMAO->liquidus.size()" ) );
        computed_quantities[q][output_index] = AO->liquidus[q];
      }
      else if ( property_names[i] == "T adiabatic" )
      {
        AssertThrow ( q < AO->Ta.size(),
          ExcMessage ( "EC_Visualization_AO::evaluate_vector_field -- q >= CMAO->Ta.size()" ) );
        computed_quantities[q][output_index] = AO->Ta[q];
      }
      else if ( property_names[i] == "T parametric" )
      {
        AssertThrow ( q < AO->Tss.size(),
          ExcMessage ( "EC_Visualization_AO::evaluate_vector_field -- q >= CMAO->Tss.size()" ) );
        computed_quantities[q][output_index] = AO->Tss[q];
      }
      else if ( property_names[i] == "F parametric" )
      {
        AssertThrow ( q < AO->Finter.size(),
          ExcMessage ( "EC_Visualization_AO::evaluate_vector_field -- q >= CMAO->Finter.size()" ) );
        computed_quantities[q][output_index] = AO->Finter[q];
      }
      else if ( property_names[i] == "T star" )
      {
        AssertThrow ( q < AO->Tstar.size(),
          ExcMessage ( "EC_Visualization_AO::evaluate_vector_field -- q >= CMAO->Tstar.size()" ) );
        computed_quantities[q][output_index] = AO->Tstar[q];
      }
    }
  }
}

template <int dim> void EC_Visualization_AO<dim>::parse_parameters ( ParameterHandler & prm )
{
  prm.enter_subsection("Postprocess");
  {
    prm.enter_subsection("Visualization");
    {
      prm.enter_subsection("Enhanced capitanio");
      {
        // Get property names and compare against variable names
        property_names = Utilities::split_string_list ( prm.get ( "List of properties" ) );
        AssertThrow(Utilities::has_unique_entries(property_names),
                    ExcMessage("The list of strings for the parameter "
                                "'Postprocess/Visualization/Capitanio Melt/List of material properties' "
                                "contains entries more than once. This is not allowed. "
                                "Please check your parameter file."));
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();
}

// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(EC_Visualization_AO,
                                                  "enhanced capitanio",
                                                  "")
    }
  }
}
