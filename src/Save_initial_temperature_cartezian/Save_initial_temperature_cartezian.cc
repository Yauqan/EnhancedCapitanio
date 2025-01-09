
#include "Save_initial_temperature_cartezian.h"
#include "../hlp/hlp.h"

#include <fstream>
#include <deal.II/numerics/vector_tools.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <sstream>

using namespace aspect;
using namespace std;
using namespace aspect::Postprocess;

template <typename T> inline bool output_now ( std::vector<T> & requested, const T & now );


template <int dim> inline void Save_initial_temperature_cartezian<dim>::make_points ()
{
  if ( volume_points.size() == 0 )
  {
    AssertThrow ( dim == 2,
      ExcMessage ( "EC_Save_initial_temperature_cartezian::make_points -- Does not yet support 3 dimensions" ) );
    
    const double dx = extents[0]/point_amounts[0];
    const double dy = extents[1]/point_amounts[1];
    for ( double x = dx/2.0; x < extents[0]; x += dx )
    {
      for ( double y = dy/2.0; y < extents[1]; y += dy )
      {
        volume_points.emplace_back ( Point<dim> ( x, y ) );
      }
    }
  }
}

template <int dim> inline void Save_initial_temperature_cartezian<dim>::fill_points ()
{
  AssertThrow ( volume_points.size() > 0,
    ExcMessage ( "EC_Save_initial_temperature_cartezian::fill_points: No made points found!" ) );

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
    AssertThrow ( n_procs > 0,
      ExcMessage (  "While trying to evaluate the solution at point " +
                    Utilities::to_string ( volume_points[p][0] ) + ", " +
                    Utilities::to_string ( volume_points[p][1] ) + "), "
                    "no processors reported that the point lies inside the "
                    "set of cells they own. Are you trying to evaluate the "
                    "solution at a point that lies outside of the domain?" ) );

    // Reduce all collected values into local Vector
    Utilities::MPI::sum ( point_values[p], this->get_mpi_communicator(), point_values[p] );

    // Normalize in cases where points are claimed by multiple processors
    if ( n_procs > 1 )
      point_values[p] /= n_procs;
  }
}

template <int dim> pair<string, string> Save_initial_temperature_cartezian<dim>::execute ( TableHandler & )
{
  if ( output_now ( output_times, this->get_time() ) == false && output_now ( output_timesteps, (int) this->get_timestep_number() ) == false )
    return {"", ""};
  
  make_points();
  fill_points();
  const int temperature_indx = 1 + dim;
  const string filename = this->get_output_directory() + "SaveInitialTemperature/temperature." + to_string ( this->get_timestep_number() );
  ofstream f;
  if ( Utilities::MPI::this_mpi_process ( this->get_mpi_communicator() ) == 0 )
  {
    f.open ( filename );
    f.setf(std::ios::fixed, std::ios::floatfield);
    f.precision(15);

    f << "# POINTS: " << point_amounts[0] << ' ' << point_amounts[1] << '\n';
  }
  AssertThrow ( dim == 2,
    ExcMessage ( "Error -- CapitanioMeltTopography::execute: Does not yet support 3 dimensions" ) );
  if ( Utilities::MPI::this_mpi_process ( this->get_mpi_communicator() ) == 0 )
    for ( int y = 0; y < point_amounts[1]; y++ )
      for ( int x = 0; x < point_amounts[0]; x++ )
        f << volume_points [ x*point_amounts[1] + y ][0] << '\t' << volume_points [ x*point_amounts[1] + y ][1] << '\t' << point_values [ x*point_amounts[1] + y ] [ temperature_indx ] << '\n';
  
  return make_pair ( string ("Writing Initial temperature values:"), filename );
}

template <int dim> void Save_initial_temperature_cartezian<dim>::declare_parameters ( ParameterHandler & prm )
{
  prm.enter_subsection ( "Postprocess" );
    prm.enter_subsection ( "Save initial temperature cartezian" );
      prm.declare_entry ( "Width extent", "1.0", Patterns::Double (0.0), "" );
      prm.declare_entry ( "Depth extent", "1.0", Patterns::Double (0.0), "" );
      prm.declare_entry ( "Width", "1", Patterns::Integer (0), "" );
      prm.declare_entry ( "Depth", "1", Patterns::Integer (0), "" );
      prm.declare_entry ( "Output times", "1e11", Patterns::List ( Patterns::Double (0.0) ), "" );
      prm.declare_entry ( "Output timesteps", to_string(numeric_limits<int>::max()), Patterns::List ( Patterns::Integer (0) ), "" );
    prm.leave_subsection ();
  prm.leave_subsection ();
}

template <int dim> void Save_initial_temperature_cartezian<dim>::parse_parameters ( ParameterHandler & prm )
{
  if ( Utilities::MPI::this_mpi_process ( this->get_mpi_communicator() ) == 0 )
    aspect::DizHlp::ensureFolderExists ( this->get_output_directory() + "SaveInitialTemperature" );
  point_amounts.resize(2,0);
  extents.resize(2,0);
  prm.enter_subsection ( "Postprocess" );
    prm.enter_subsection ( "Save initial temperature cartezian" );
      point_amounts[0]      = prm.get_integer ( "Width" );
      point_amounts[1]      = prm.get_integer ( "Depth" );
      extents[0]            = prm.get_double ( "Width extent" );
      extents[1]            = prm.get_double ( "Depth extent" );
      output_times          = Utilities::string_to_double ( Utilities::split_string_list ( prm.get ( "Output times" ) ) );
      if ( this->convert_output_to_years() )
        for ( double & ot: output_times )
          ot *= year_in_seconds;
      output_timesteps      = Utilities::string_to_int ( Utilities::split_string_list ( prm.get ( "Output timesteps" ) ) );
      outputted_times       = vector<bool> ( output_times.size(), false );
      outputted_timesteps   = vector<bool> ( output_timesteps.size(), false );
    prm.leave_subsection ();
  prm.leave_subsection ();
}

template <int dim> template <class Archive> void Save_initial_temperature_cartezian<dim>::serialize ( Archive & ar, const unsigned int )
{
  ar & output_times
  & output_timesteps
  ;
}

template <int dim> void Save_initial_temperature_cartezian<dim>::save ( map<string, string> & status_strings ) const
{
  std::ostringstream os;
  {
    aspect::oarchive oa (os);
    oa << (*this);
  }
  status_strings["Save_initial_temperature_cartezian"] = os.str();
}

template <int dim> void Save_initial_temperature_cartezian<dim>::load ( const map<string, string> & status_strings )
{
  if (status_strings.find("Save_initial_temperature_cartezian") != status_strings.end())
  {
    std::istringstream is (status_strings.find("Save_initial_temperature_cartezian")->second);
    aspect::iarchive ia (is);
    ia >> (*this);
  }
}

namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR ( Save_initial_temperature_cartezian,
                                    "save initial temperature cartezian",
                                    "" )
  }
}





















template <typename T> inline bool output_now ( std::vector<T> & requested, const T & now )
{
  for ( long unsigned int k = 0; k < requested.size(); ++k )
    if ( requested[k] <= now )
    {
      requested.erase ( requested.begin() + k );
      return true;
    }
  return false;
}
