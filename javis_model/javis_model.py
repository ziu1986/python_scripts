class JavisModel:
    def __init__(self, species, **karg):
        self.name = species
        self.temp_min = karg.pop("Tmin")
        self.temp_max = karg.pop("Tmax")
        self.temp_opt = karg.pop("Topt")
        self.f_min = karg.pop("fmin")
        self.vpd_min = karg.pop("Dmin")
        self.vpd_max = karg.pop("Dmax")
        self.alpha_light = karg.pop("alpha")
        self.gmax = karg.pop("gmax") # mmol O3 m-2 PLS s-1
        self.sgs = karg.pop('start_gs', 0)
        self.egs = karg.pop('end_gs', 366)


    def f_temp(self, temperature):
        '''
        Compute Javis f_temp function.
        Parameters
        ----------
        temperature : float
        The 2m temperature in deg C
        Keyword arguments
        -----------------
        Tmin : float
        Minimum temperature for stomatal opening in deg C
        Tmax : float
        Maximum temperature for stomatal opening in deg C
        Topt : float
        Optimal temperature for stomatal opening in deg C
        Returns
        -------
        f_temp : float, range 0,1
        '''
        import numpy as np
       
        TMIN = self.temp_min
        TMAX = self.temp_max
        TOPT = self.temp_opt
        FMIN = self.f_min

        # Catch tamperature > TMAX
        try:
            temperature[np.where(temperature>TMAX)] = TMAX
        except TypeError:
            if temperature > TMAX:
                temperature = TMAX
                
        # Compute f_temp function
        beta = (TMAX-TOPT)/(TOPT-TMIN)
        f_temp = (temperature-TMIN)/(TOPT-TMIN)*((TMAX-temperature)/(TMAX-TOPT))**beta

        # Catch invalid values
        try:
            f_temp[np.where(f_temp<FMIN)] = FMIN
        except TypeError:
            if f_temp < FMIN:
                f_temp = FMIN
        # Catch values larger 1
        try:
            f_temp[np.where(f_temp>1)] = 1
        except TypeError:
            if f_temp > 1:
                f_temp = 1

        return(f_temp) 

    def f_vpd(self, vpd):
        '''
        Compute Javis f_vpd function.
        Parameters
        ----------
        VPD : float
        Water vapor pressure deficit in kPa
        Keyword arguments
        -----------------
        fmin : float
        Minimum value for stomatal opening
        Dmin : float
        Maximum VPD for stomatal opening in KPa
        Dmax : float
        Maximum VPD for stomatal opening in kPa
        Returns
        -------
        f_vpd : float, range 0,1
        
        '''
        import numpy as np
        
        FMIN = self.f_min
        DMIN = self.vpd_min
        DMAX = self.vpd_max

        f_vpd = FMIN + (1-FMIN) * (DMIN-vpd)/(DMIN-DMAX)

        # Catch negative values
        try:
            f_vpd[np.where(f_vpd<FMIN)] = FMIN
        except TypeError:
            if f_vpd < FMIN:
                f_vpd = FMIN
        # Catch values larger 1
        try:
            f_vpd[np.where(f_vpd>1)] = 1
        except TypeError:
            if f_vpd > 1: 
                f_vpd = 1

        return(f_vpd)


    def f_light(self, ppfd):
        '''
        Compute Javis f_light function.
        Parameters
        ----------
        ppfd : float
        Photoactive photonflaux density in W/m^2/s
        Keyword arguments
        -----------------
        alpha : float
        slope parameter in m^2*s/W
        Returns
        -------
        f_light : float, range 0,1
        
        '''

        import numpy as np
        
        ALPHA = self.alpha_light

        f_light = 1 - np.exp(-ALPHA*ppfd)
        
        # Catch negative values
        try:
            f_light[np.where(f_light<0)] = 0
        except TypeError:
            if f_light < 0:
                f_light = 0
        # Catch values larger 1
        try:
            f_light[np.where(f_light>1)] = 1
        except TypeError:
            if f_light >1:
                f_light = 1

        return(f_light)
    
    def f_sw(self, soil_inst):
        '''
        Compute the soilwater stress by plant available water (PAW) method.
        This will need a soil model.
        Parameters
        ----------
        FC : float
            Soil caracteristics
        PAW : float
            plant available water
        Returns
        -------
        f_sw : float

        '''
        return(0)

    def f_phen(self, date):
        '''
        Compute the phenology of the species.
        Parameters
        ----------
        Date : datetime64
        Arguments
        -----------------
        start_gs : int
            Start of the growing season in day of thew year.
        end_gs : int
            End of the growing season in day of the year.
        Returns
        -------
        f_phen : float
            Stage of the phenology.
        '''
        import numpy as np
        f_phen = np.zeros(date.size)
        f_phen[np.where((date.dayofyear>=start_gs) & (data.dayofyear<=end_gs))] = 1

        return(f_phen)
        
        
    def stomatal_conductance(self, temperature, vpd, ppfd, date):
        '''
        Compute Javis stomatal conductance.
        Parameters
        ----------
        f_phen : float
            
        f_light : float
        
        f_min : float
        
        f_temp : float
        
        f_vpd : float
        
        f_sw : float
            
        Arguments
        ---------
        gmax : float
        Maximum stomatal conductance
        Returns
        -------
        gsto : float
            Stomatal conductance

        '''
        # Compute f-functions
        v_temp = self.f_temp(temperature)
        v_light = self.f_light(ppfd)
        v_vpd = self.f_vpd(vpd)
        v_min = self.f_min
        v_phen = self.f_phen(date)
        
        gsto = self.gmax * v_phen * v_light * np.maximum(v_min, v_temp * v_vpd * v_sw)

        return(gsto)





    
