
#include "EC-Melt_topography_postprocessor.h"
#include "../hlp/hlp.h"

#include <fstream>
#include <deal.II/numerics/vector_tools.h>

using namespace aspect;
using namespace std;
using namespace aspect::Postprocess;

template <int dim> inline void EC_Melt_topography<dim>::make_points ()
{
  if ( volume_points.size() == 0 )
  {
    AssertThrow ( dim == 2,
      ExcMessage ( "Does not yet support 3 dimensions" ) );
    const double dx = extents[0]/point_amounts[0];
    dy = extents[1]/point_amounts[1];
    for ( double x = dx/2.0; x < extents[0]; x += dx )
      for ( double y = dy/2.0; y < extents[1]; y += dy )
        volume_points.emplace_back ( Point<dim> ( x, y ) );
  }
  // if ( surface_points.size() == 0 )
  // {
  //   AssertThrow ( dim == 2,
  //     ExcMessage ( "Does not yet support 3 dimensions" ) );
  //   const double dx = extents[0]/point_amounts[0];
  //   for ( double x = dx/2.0; x < extents[0]; x += dx )
  //     surface_points.emplace_back ( Point<dim> ( x, 0.0 ) );
  // }
}

template <int dim> inline void EC_Melt_topography<dim>::fill_points ()
{
  AssertThrow ( volume_points.size() > 0,
                ExcMessage ( "No made points found!" ) );

  point_values = vector<Vector<double>> ( volume_points.size(), Vector<double> ( this->introspection().n_components ) );
  for ( unsigned int p = 0; p < volume_points.size(); ++p )
  {
    bool point_found = false;
    try
    {
      VectorTools::point_value ( this->get_mapping(),
                                 this->get_dof_handler(),
                                 this->get_solution(),
                                 volume_points[p],
                                 point_values[p] );
      point_found = true;
    }
    catch ( const VectorTools::ExcPointNotAvailableHere & )
    {;}

    // ensure that at least one processor found things
    const int n_procs = Utilities::MPI::sum ( point_found ? 1 : 0, this->get_mpi_communicator() );
    AssertThrow ( n_procs > 0, ExcMessage ( "While trying to evaluate the solution at point " +
                                            Utilities::to_string ( volume_points[p][0] ) + ", " +
                                            Utilities::to_string ( volume_points[p][1] ) +
                                            + "), " +
                                            "no processors reported that the point lies inside the " +
                                            "set of cells they own. Are you trying to evaluate the " +
                                            "solution at a point that lies outside of the domain?"
                                          ) );

    // Reduce all collected values into local Vector
    Utilities::MPI::sum ( point_values[p], this->get_mpi_communicator(), point_values[p] );

    // Normalize in cases where points are claimed by multiple processors
    if ( n_procs > 1 )
      point_values[p] /= n_procs;
  }

}

template <int dim> pair<string, string> EC_Melt_topography<dim>::execute ( TableHandler & )
{
  if ( std::isnan ( last_output_time ) )
  {
    last_output_time = this->get_time() - output_interval;
    last_output_timestep = this->get_timestep_number();
  }
  if ( ( this->get_time() < last_output_time + output_interval )
      && ( this->get_timestep_number() < last_output_timestep + maximum_timesteps_between_outputs )
      && ( this->get_timestep_number() != 0 ) )
    return {"", ""};
  
  make_points();
  fill_points();
  const int capitanio_melt_indx = this->introspection().compositional_index_for_name("mantle_depletion") + 2 + dim;
  const string filename = this->get_output_directory() + "Topography/melt_topography." + Utilities::int_to_string(output_index, 4);
  ofstream f;
  if ( Utilities::MPI::this_mpi_process ( this->get_mpi_communicator() ) == 0 )
    f.open ( filename );
  
  AssertThrow ( dim == 2, 
    ExcMessage ( "Does not yet support 3 dimensions" ) );
  for ( int x = 0; x < point_amounts[0]; x++ )
  {
    double capitanio_melt = 0.0;
    for ( int y = 0; y < point_amounts[1]; y++ )
      capitanio_melt += point_values [ x*point_amounts[1] + y ] [ capitanio_melt_indx ];
    const double h_C = Utilities::MPI::sum ( capitanio_melt, this->get_mpi_communicator() ) * dy;
    if ( Utilities::MPI::this_mpi_process ( this->get_mpi_communicator() ) == 0 )
    {
      f << volume_points [ x*point_amounts[1] ][0] << '\t'
      << h_C * ( 1.0 - rho_melt / rho_lit )
      << '\n';
    }
  }
  
  output_index ++;
  set_last_output_time ( this->get_time() );
  last_output_timestep = this->get_timestep_number();
  return make_pair ( string ("Writing Enhanced capitanio melt topography values:"), filename );
}

template <int dim> void EC_Melt_topography<dim>::declare_parameters ( ParameterHandler & prm )
{
  prm.enter_subsection ( "Postprocess" );
    prm.enter_subsection ( "Enhanced capitanio melt topography" );
      prm.declare_entry ( "Width extent", "1.0", Patterns::Double (0.0), "" );
      prm.declare_entry ( "Depth extent", "1.0", Patterns::Double (0.0), "" );
      prm.declare_entry ( "Width", "1", Patterns::Integer (0), "" );
      prm.declare_entry ( "Depth", "1", Patterns::Integer (0), "" );
      prm.enter_subsection ( "Physical parameters" );
        prm.declare_entry ( "Melt density", "2600", Patterns::Double (0.0), "Units: \\si{\\kilogram\\meter^3}" );
        prm.declare_entry ( "Lithosphere density", "3240", Patterns::Double (0.0), "Units: \\si{\\kilogram\\meter^3}" );
      prm.leave_subsection ();
    prm.leave_subsection ();
  prm.leave_subsection ();
}

template <int dim> void EC_Melt_topography<dim>::parse_parameters ( ParameterHandler & prm )
{
  if ( Utilities::MPI::this_mpi_process ( this->get_mpi_communicator() ) == 0 )
    aspect::DizHlp::ensureFolderExists ( this->get_output_directory() + "Topography" );
  output_index = 0;

  point_amounts.resize(2,0);
  extents.resize(2,0);
  prm.enter_subsection ( "Postprocess" );
    prm.enter_subsection ( "Enhanced capitanio melt topography" );
      point_amounts[0]            = prm.get_integer ( "Width" );
      point_amounts[1]            = prm.get_integer ( "Depth" );
      extents[0]                  = prm.get_double ( "Width extent" );
      extents[1]                  = prm.get_double ( "Depth extent" );
      prm.enter_subsection ( "Physical parameters" );
        rho_melt                  = prm.get_double ( "Melt density" );
        rho_lit                   = prm.get_double ( "Lithosphere density" );
      prm.leave_subsection ();
    prm.leave_subsection ();
    prm.enter_subsection ( "Visualization" );
      output_interval = prm.get_double ( "Time between graphical output" );
      if ( this->convert_output_to_years() )
        output_interval *= year_in_seconds;
      maximum_timesteps_between_outputs = prm.get_integer ( "Time steps between graphical output" );
      if (output_interval > 0.0)
        AssertThrow ( this->get_parameters().run_postprocessors_on_nonlinear_iterations == false,
                      ExcMessage ( "Postprocessing nonlinear iterations is only supported if every time "
                                   "step is visualized, or in other words, if the 'Time between graphical "
                                   "output' in the Visualization postprocessor is set to zero." ) );
    prm.leave_subsection ();
  prm.leave_subsection ();
}

template <int dim> void EC_Melt_topography<dim>::set_last_output_time ( const double current_time )
{
  if ( output_interval > 0 )
  {
    const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
    last_output_time = last_output_time + std::floor((current_time-last_output_time)/output_interval*magic) * output_interval/magic;
  }
}

template <int dim> template <class Archive> void EC_Melt_topography<dim>::serialize ( Archive & ar, const unsigned int )
{
  ar & last_output_time
  & last_output_timestep
  & output_index
  ;
}

template <int dim> void EC_Melt_topography<dim>::save ( map<string, string> & status_strings ) const
{
  std::ostringstream os;
  {
    aspect::oarchive oa (os);
    oa << (*this);
  }
  status_strings["CapitanioMeltTopography"] = os.str();
}

template <int dim> void EC_Melt_topography<dim>::load ( const map<string, string> & status_strings )
{
  if (status_strings.find("CapitanioMeltTopography") != status_strings.end())
  {
    std::istringstream is (status_strings.find("CapitanioMeltTopography")->second);
    aspect::iarchive ia (is);
    ia >> (*this);
  }
}

namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR ( EC_Melt_topography,
                                    "enhanced capitanio melt topography",
                                    "" )
  }
}
