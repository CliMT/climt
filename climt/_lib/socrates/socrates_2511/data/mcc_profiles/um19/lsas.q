netcdf lsas.q{                                                                  

dimensions:
    lat            =   1;
    lon            =   1;
    plev           =  19;


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

    float q(plev,lon,lat);                                                     
             q:units = "None";                                                
             q:title = "MMR OF WATER VAPOUR";                                                

data:                                                                           
              lat =   .000000E+00;
              lon =   .000000E+00;
             plev =   .530250E+03,  .151500E+04,  .303000E+04,  .580750E+04,
                      .101000E+05,  .151500E+05,  .202000E+05,  .252500E+05,
                      .303000E+05,  .358550E+05,  .426725E+05,  .510050E+05,
                      .606000E+05,  .707000E+05,  .800425E+05,  .878699E+05,
                      .939804E+05,  .984750E+05,  .100697E+06;
                q =   .399748E-05,  .399885E-05,  .400499E-05,  .399914E-05,
                      .403151E-05,  .515147E-05,  .148531E-04,  .341285E-04,
                      .804108E-04,  .245426E-03,  .560977E-03,  .108482E-02,
                      .193221E-02,  .307083E-02,  .428539E-02,  .517420E-02,
                      .610843E-02,  .693940E-02,  .739454E-02;

}                                                                               
