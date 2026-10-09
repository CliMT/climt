netcdf lsaw.t{                                                                  

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
             plev =   .531825E+03,  .151950E+04,  .303900E+04,  .582475E+04,
                      .101300E+05,  .151950E+05,  .202600E+05,  .253250E+05,
                      .303900E+05,  .359615E+05,  .427993E+05,  .511565E+05,
                      .607800E+05,  .709100E+05,  .802802E+05,  .881309E+05,
                      .942596E+05,  .987675E+05,  .100996E+06;
                t =   .215629E+03,  .212847E+03,  .213643E+03,  .215894E+03,
                      .217026E+03,  .217200E+03,  .217200E+03,  .217932E+03,
                      .221402E+03,  .228107E+03,  .236307E+03,  .244609E+03,
                      .251201E+03,  .255221E+03,  .257818E+03,  .258547E+03,
                      .257832E+03,  .257313E+03,  .257123E+03;

}                                                                               
