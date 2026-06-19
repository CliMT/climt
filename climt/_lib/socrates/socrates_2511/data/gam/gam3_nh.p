netcdf gam3_nh.p{                                                               

dimensions:
    lat            =   1;
    lon            =   1;
    plev           =  18;


variables:
    float lat(lat);                                                           
             lat:units = "degree";                                            
             lat:title = "LATITUDE";                                         
    float lon(lon);                                                           
             lon:units = "degree";                                            
             lon:title = "LONGITUDE";                                        
    float plev(plev);                                                         
             plev:units = "Pa";                                               
             plev:title = "PRESSURE";                                        

    float p(plev,lon,lat);                                                     
             p:units = "Pa";                                                  
             p:title = "PRESSURE";                                                           

data:                                                                           
              lat =   .000000E+00;
              lon =   .000000E+00;
             plev =   .100000E+03,  .300000E+03,  .100000E+04,  .300000E+04,
                      .500000E+04,  .700000E+04,  .100000E+05,  .150000E+05,
                      .200000E+05,  .202088E+05,  .204119E+05,  .250000E+05,
                      .300000E+05,  .400000E+05,  .500000E+05,  .700000E+05,
                      .850000E+05,  .944467E+05;
                p =   .100000E+03,  .300000E+03,  .100000E+04,  .300000E+04,
                      .500000E+04,  .700000E+04,  .100000E+05,  .150000E+05,
                      .200000E+05,  .202088E+05,  .204119E+05,  .250000E+05,
                      .300000E+05,  .400000E+05,  .500000E+05,  .700000E+05,
                      .850000E+05,  .944467E+05;

}                                                                               
