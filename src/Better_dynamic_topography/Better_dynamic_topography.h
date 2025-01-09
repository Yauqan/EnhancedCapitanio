
#ifndef _aspect_postprocess_better_surface_topography_h
#define _aspect_postprocess_better_surface_topography_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    class BetterDynamicTopography : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        std::pair<std::string,std::string> execute ( TableHandler & statistics ) override;
        const LinearAlgebra::BlockVector & topography_vector() const;
        const Vector<float> & cellwise_topography() const;

        std::list<std::string> required_other_postprocessors() const override;
        void parse_parameters ( ParameterHandler & prm ) override;
        static void declare_parameters ( ParameterHandler & prm );

        template <class Archive> void serialize ( Archive & ar, const unsigned int version );
        void save ( std::map<std::string, std::string> & status_strings ) const override;
        void load ( const std::map<std::string, std::string> & status_strings ) override;
      private:
        void output_to_file ( bool upper, std::vector<std::pair<Point<dim>, double>> & values );

        void set_last_output_time ( const double current_time );

        LinearAlgebra::BlockVector topo_vector;
        Vector<float> visualization_values;
        double density_above;
        double density_below;
        bool output_surface;
        bool output_bottom;

        unsigned int maximum_timesteps_between_outputs;
        double output_interval;
        unsigned int last_output_timestep = numbers::invalid_unsigned_int;
        double last_output_time = std::numeric_limits<double>::quiet_NaN();
        int output_index;
    };
  }
}


#endif
