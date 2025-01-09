

#ifndef _aspect_enhanced_capitanio_material_model_h
#define _aspect_enhanced_capitanio_material_model_h

#include <aspect/material_model/interface.h>
#include <aspect/material_model/rheology/drucker_prager.h>
#include <aspect/material_model/equation_of_state/linearized_incompressible.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator_access.h>

#include <vector>
#include <string>

#include "../EC-Additional_outputs/EC-Additional_outputs.h"

namespace aspect
{
  namespace MaterialModel
  {
      template <int dim> class EC_Material_model : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:
          void evaluate ( const MaterialModel::MaterialModelInputs<dim> & in, MaterialModel::MaterialModelOutputs<dim> & out ) const override;
          bool is_compressible () const override;
          static void declare_parameters ( ParameterHandler & prm );
          void create_additional_named_outputs ( MaterialModelOutputs<dim> & out ) const override;
          void parse_parameters ( ParameterHandler & prm ) override;

        private:
          double thermal_conductivity;
          EquationOfState::LinearizedIncompressible<dim> equation_of_state;

          struct melt_parameters
          {
            bool calculate_depletion;
            bool use_depth_regeneration;
            double regeneration_depth;
          } par_melt;

          struct viscosity_parameters
          {
            bool use_druckerprager;
            bool use_melt_depletion;
            bool use_diffusion;
            double eta_0;
            struct diffusion_struct
            {
              double A;
              double E_a;
            } diffusion;
            struct limits_struct
            {
              double eta_min;
              double eta_max;
            } limits;
            struct depletion_struct
            {
              double k;
            } depletion;
            struct druckerprager_struct
            {
              double C;
              double mu;
              double em;
              bool use_melt_depletion;
              double delta_C_deplt;
              double delta_mu_deplt;
            } druckerprager;
          } par_viscosity;

          struct density_parameters
          {
            bool use_melt_depletion;
            struct depletion_struct
            {
              double delta_rho;
            } depletion;
          } par_density;

          const double R = 8.314;     // Molarni plynova konstanta
      };
  }
}

#endif
