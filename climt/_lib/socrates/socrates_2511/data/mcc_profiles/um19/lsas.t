netcdf lsas.t{                                                                  

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

    float t(plev,lon,lat);                                                     
             t:units = "K";                                                   
             t:title = "TEMPERATURE";                                                        

data:                                                                           
              lat =   .000000E+00;
              lon =   .000000E+00;
             plev =   .530250E+03,  .151500E+04,  .303000E+04,  .580750E+04,
                      .101000E+05,  .151500E+05,  .202000E+05,  .252500E+05,
                      .303000E+05,  .358550E+05,  .426725E+05,  .510050E+05,
                      .606000E+05,  .707000E+05,  .800425E+05,  .878699E+05,
                      .939804E+05,  .984750E+05,  .100697E+06;
                t =   .238554E+03,  .230013E+03,  .225896E+03,  .225000E+03,
                      .225000E+03,  .225000E+03,  .225000E+03,  .227906E+03,
                      .234790E+03,  .243274E+03,  .252323E+03,  .261143E+03,
                      .268385E+03,  .273990E+03,  .278794E+03,  .282531E+03,
                      .284969E+03,  .286408E+03,  .286937E+03;

}                                                                               
