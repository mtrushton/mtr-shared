def calc_chisq(flux_data, flux_data_error, flux_model):
    """
    Scale FLUX_MODEL to FLUX_DATA and calculate chisq.
    
    Inputs
    -----------
    flux_data (array): Observed flux 
    flux_data_error(array): Error in the observed data
    flux_model(array) : Model flux
        
    Outputs
    --------
    chisq (float): chi-square statistic
    flux_scaled_model(array) : Model flux scaled to FLUX_DATA
    """
    # Calculate scaling parameters
    sx = sum(flux_model / flux_data_error**2)
    sy = sum(flux_data / flux_data_error**2)
    schi = sum(1 / flux_data_error**2)

    tchi = (flux_model - sx / schi) / flux_data_error

    stt = sum(tchi**2)
    ty = sum(tchi * flux_data / flux_data_error)

    bchi = ty / stt
    achi = ((sy - (sx * bchi)) / schi)

    # Scale the model
    flux_scaled_model = achi + bchi * flux_model
    
    # Calculate chisq
    chisq = sum(((flux_data - flux_scaled_model) / flux_data_error)**2)
    
    return chisq, flux_scaled_model
