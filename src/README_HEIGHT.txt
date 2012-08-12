"""
        # Estimate canopy height from stem C, from JULES
        a_ws = 10.0   # woody biomass as a multiply of live stem biomass: 1.0 = grass, 10 = shrub
        eta_sl = 0.01 # live stemwood coefficient
        b_wl = 1.667  # allometric exponent relating woody biomass to the lai 1.667
        a_wl = 0.65   # allometric coeff relating woody biomass to lai: 0.005 = grass, 0.1 = shrub
        tonnes_per_ha_2_g_m2 = 0.01
        g_to_kg = 0.001
        stem = self.state.stem / tonnes_per_ha_2_g_m2 * g_to_kg
        height = stem / (a_ws * eta_sl) * (a_wl / stem)**(1.0 / b_wl)
        print height
        """