
#ifndef _aspect_postprocess_save_initial_temperature_cartezian_h
#define _aspect_postprocess_save_initial_temperature_cartezian_h

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
    template <int dim> class Save_initial_temperature_cartezian : public Interface<dim>, public SimulatorAccess<dim>
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

        std::vector<int> point_amounts;
        std::vector<double> extents;
        std::vector<int> output_timesteps;
        std::vector<double> output_times;
        std::vector<bool> outputted_times;
        std::vector<bool> outputted_timesteps;

        std::vector<Point<dim>> volume_points;
        std::vector<Vector<double>> point_values;
    };
  }
}

#endif
