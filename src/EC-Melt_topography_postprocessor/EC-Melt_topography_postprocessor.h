
#ifndef _aspect_postprocess_enhanced_capitanio_melt_topography_h
#define _aspect_postprocess_enhanced_capitanio_melt_topography_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/data_out_base.h>

#include <utility>
#include <vector>
#include <string>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim> class EC_Melt_topography : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        std::pair<std::string, std::string> execute ( TableHandler & statistics ) override;
        static void declare_parameters ( ParameterHandler & prm );
        void parse_parameters ( ParameterHandler & prm ) override;
        template <class Archive> void serialize ( Archive & ar, const unsigned int version );
        void save ( std::map<std::string, std::string> & status_strings ) const override;
        void load ( const std::map<std::string, std::string> & status_strings ) override;
      private:
        inline void make_points ();
        inline void fill_points ();

        void set_last_output_time ( const double current_time );

        std::vector<int> point_amounts;
        std::vector<double> extents;
        unsigned int maximum_timesteps_between_outputs;
        double output_interval;
        unsigned int last_output_timestep = numbers::invalid_unsigned_int;
        double last_output_time = std::numeric_limits<double>::quiet_NaN();
        int output_index;
        
        double dy;

        // double delta_rho;
        double rho_melt;
        // double h_lit;
        double rho_lit;
        // double eta_lit;

        std::vector<Point<dim>> volume_points;
        // std::vector<Point<dim>> surface_points;
        std::vector<Vector<double>> point_values;
        // std::vector<Vector<double>> surface_point_values;
        // std::vector<Vector<double>> surface_point_gradient;
    };
  }
}

#endif
