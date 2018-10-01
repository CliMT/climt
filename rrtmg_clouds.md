# RRMTG Clouds with McICA

## Cloud properties

There are three options for the RRTMG `inflag`, as given in the `climt` dictionary:

```python
rrtmg_cloud_props_dict = {
    'direct_input': 0,
    'single_cloud_type': 1,
    'liquid_and_ice_clouds': 2
}
```

For McICA, we cannot use `single_cloud_type`, but can choose between `direct_input` and `liquid_and_ice_clouds`.
If, we choose `direct_input`, the cloud ice and water droplet particle sizes do not do anything, but the input cloud optical depth matters (for both the long and shortwave).
Similarly, for the shortwave, the `single_scattering_albedo_due_to_cloud`, `cloud_asymmetry_parameter` and `cloud_forward_scattering_fraction` matter.
Here the RRTMG `iceflag` and `liqflag` do not matter.

On the other hand, if we choose `liquid_and_ice_clouds`, the optical depth is calculated from the cloud ice and water droplet particle sizes (`cloud_ice_particle_size` and `cloud_water_droplet_radius`), as well as the cloud ice and water paths (`mass_content_of_cloud_ice_in_atmosphere_layer`, `mass_content_of_cloud_liquid_water_in_atmosphere_layer`).
Likewise, the single scattering albedo and the cloud asymmetry parameter (both for the shortwave) are calculated from the cloud ice and water properties in the case of `liquid_and_ice_clouds`.
As such, any input values for `longwave_optical_thickness_due_to_cloud`, `shortwave_optical_thickness_due_to_cloud`, `single_scattering_albedo_due_to_cloud`, `cloud_asymmetry_parameter` and `cloud_forward_scattering_fraction` are irrelevant.
The calculation of the above variables depends on the cloud ice and water properties; `iceflag` and `liqflag` in RRTMG.


## Calculation of cloud properties

### Longwave

For the longwave, the only cloud property of interest for calculating radiative fluxes in RRTMG, is the optical depth. This is calculated as follows.

```fortran
do lay = 1, nlayers
  do ig = 1, ngptlw
    taucmc(ig,lay) = ciwpmc(ig,lay) * abscoice(ig) + &
                     clwpmc(ig,lay) * abscoliq(ig)
  end do
end do
```

Values of cloud optical depth `taucmc` are calculated for each model layer (pressure), `lay`, and each g-interval, `ig`.
The cloud ice and liquid absorption coefficients (`abscoice` and `abscoliq`) are multiplied by the cloud ice and liquid water paths (`ciwpmc` and `clwpmc`) respectively, to give the ice cloud optical depth and the liquid water cloud optical depth.
The cloud ice and liquid water paths are input by the user, in `climt` as `mass_content_of_cloud_ice_in_atmosphere_layer` and `mass_content_of_cloud_liquid_water_in_atmosphere_layer` respectively.
The cloud ice and liquid absorption coefficients are calculated based on the ice and liquid water particle sizes (specified by the user), and this calculation depends on the choice of `iceflag` and `liqflag`.

### Shortwave

For the shortwave, there are three cloud properties, which affect the radiative flux calculation in RRTMG, namely the optical depth, the single scattering albedo and the asymmetry parameter.

The shortwave optical depth is calculated as:
```fortran
taucmc(ig,lay) = tauice + tauliq
```
with
```fortran
tauice = (1 - forwice(ig) + ssacoice(ig)) * ciwpmc(ig,lay) * extcoice(ig)
tauliq = (1 - forwliq(ig) + ssacoliq(ig)) * clwpmc(ig,lay) * extcoliq(ig)
```

The single scattering albedo is calculated as:
```fortran
ssacmc(ig,lay) = (scatice + scatliq) / taucmc(ig,lay)
```
with
```fortran
scatice = ssacoice(ig) * (1._rb - forwice(ig)) / (1._rb - forwice(ig) * ssacoice(ig)) * tauice
scatliq = ssacoliq(ig) * (1._rb - forwliq(ig)) / (1._rb - forwliq(ig) * ssacoliq(ig)) * tauliq
```

The asymmetry parameter is given by:
```fortran
asmcmc(ig,lay) = 1.0_rb / (scatliq + scatice) * ( &
                   scatliq * (gliq(ig) - forwliq(ig)) / (1.0_rb - forwliq(ig)) + &
                   scatice * (gice(ig) - forwice(ig)) / (1.0_rb - forwice(ig)) )
```

Values of optical depth, single scattering albedo and asymmetry parameter are calculated for each model layer (pressure), `lay`, and each g-interval, `ig`. The cloud ice and liquid water paths (`ciwpmc` and `clwpmc`) are input by the user.
The other parameters (`extcoice`, `extcoliq`, `ssacoice`, `ssacoliq`, `gice`, `gliq`, `forwice`, `forwliq`) are calculated based on the ice and liquid water particle sizes and this calculation depends on the choice of `iceflag` and `liqflag`.

## Cloud ice properties

There are four options for the RRTMG `iceflag`. These are given in the `climt` dictionary:

```python
rrtmg_cloud_ice_props_dict = {
    'ebert_curry_one': 0,
    'ebert_curry_two': 1,
    'key_streamer_manual': 2,
    'fu': 3
}
```

### ebert_curry_one

For the longwave, `ebert_curry_one` gives an absorption coefficient of
```fortran
abscoice(ig) = 0.005_rb + 1.0_rb/radice
```
Here, `radice` is the ice particle size and the absorption coefficient is the same for all wavebands.

`ebert_curry_one` should not be used for the shortwave component with McICA.

### ebert_curry_two

`ebert_curry_two` is the default choice for cloud optical properties in `climt`.
In this case, the longwave absorption coefficient is calculated as follows.
```fortran
do ig = 1, ngptlw
  ib = icb(ngb(ig))
  abscoice(ig) = absice1(1,ib) + absice1(2,ib)/radice
```
The absorption coefficient `abscoice` is a function of g-interval `ig` and is made up of two contributions.
The first of these `absice1(1, ib)` comes from a look up table and is given in [m<sup>2</sup> g].
`ib` provides an index for the look up table, based on the waveband of the g-interval.
`absice1(2,ib)` also comes from a look up table and is given in [microns m<sup>2</sup> g]. It is divided by `radice`, the cloud ice particle size, providing an ice particle size dependence of the absorption coefficient.
Although the syntax does not emphasise it, the absorption coefficient may also depend on model layer (pressure), as `radice` can have model layer dependence.
`radice` comes from the input property labeled `cloud_ice_particle_size` in `climt`.
The ice particle size dependent term is more important than the independent term (`absice1(2, ib)/radice` > `absice(ig)`) at all wavebands for ice particle sizes less than 88 microns.
Using `ebert_curry_two`, the ice particle size must be in the range [13, 130] microns, and even for larger particle sizes (> 88), the ice particle size dependent term is more important than the independent term for four of the five wavebands.

For the shortwave, the parameters (required for the optical depth, single scattering albedo and asymmetry) are calculated as follows.
```fortran
extcoice(ig) = abari(icx) + bbari(icx)/radice
ssacoice(ig) = 1._rb - cbari(icx) - dbari(icx) * radice
gice(ig) = ebari(icx) + fbari(icx) * radice
forwice(ig) = gice(ig)*gice(ig)
```
`abari`, `bbari`, `cbari`, `dbari`, `ebari` and `fbari` are all look up tables, containing five values, which correspond to five different wavebands.
The choice of waveband is indicated by `icx`.
The particle size dependence comes from `radice`, so each parameter consists of both a size independent and a size dependent contribution.


### key_streamer_manual

The longwave absorption coefficient is calculated as follows:

```fortran
abscoice(ig) = absice2(index,ib) + fint * (absice2(index+1,ib) - (absice2(index,ib)))
```

In this case, `index` provides a particle size dependent index for the look up table and `ib` provides a waveband dependent index.
There are 42 values for each of the 16 different wavebands.
The ice particle size must be in the range [5, 131] microns and `index` maps this size onto an integer value in the range [1, 42].
The mapping works such that particles sizes of 5, 6, and 7 microns are all mapped to 1, particles sizes of 8, 9 and 10 microns are all mapped to 2 and so on.
The second term contributing to `abscoice` is a correction term.
`fint` is between zero and one, so this term differentiates the absorption coefficient between similarly sized particles.

For the shortwave, the parameters are calculated as follows.

```fortran
extcoice(ig) = extice2(index,ib) + fint * (extice2(index+1,ib) -  extice2(index,ib))
ssacoice(ig) = ssaice2(index,ib) + fint * (ssaice2(index+1,ib) -  ssaice2(index,ib))
gice(ig) = asyice2(index,ib) + fint * (asyice2(index+1,ib) -  asyice2(index,ib))
forwice(ig) = gice(ig)*gice(ig)
```

As for the longwave calculation, `index` provides a particle size dependent index for the look up table and `ib` provides a waveband dependent index.
There are 42 values for each of the 14 different wavebands.
The mapping of particle size onto `index` is the same as for the longwave calculation.
Here again, the second term in the equations for `extcoice`, `ssacoice` and `gice` is a correction term, with `fint` between zero and one.


### fu

The longwave absorption coefficient is calculated in the same way as in `key_streamer_manual`, but the look up tables differ.
Comments in the RRTMG code state that the look up tables for `fu` are for a hexagonal ice particle parameterisation, whereas those for `key_streamer_manual` are for a spherical ice particle parameterisation.
The look up tables for `fu` are slightly larger, and the range of allowed values for the ice particle size is corresponding larger ([5, 140] microns).

For the shortwave, the coefficients `extcoice`, `ssacoice` are calculated as in `key_streamer_manual`, but from different look up tables (allowing for particle sizes in the range [5, 140] microns).
However, the calculation of `forwice` differs and it no longer depends on `gice`.

```fortran
fdelta(ig) = fdlice3(index,ib) + fint * (fdlice3(index+1,ib) - fdlice3(index,ib))
forwice(ig) = fdelta(ig) + 0.5_rb / ssacoice(ig)
```


## Cloud liquid properties

There are two options for the RRTMG `liqflag`. These are given in the `climt` dictionary:

```python
rrtmg_cloud_liquid_props_dict = {
    'radius_independent_absorption': 0,
    'radius_dependent_absorption': 1
}
```

For `radius_independent_absorption`, the longwave absorption coefficient is 0.0903614 for all wavebands.
This option should not be used for the shortwave.

For `radius_dependent_absorption`, the longwave absorption coefficient is calculated as:

```fortran
abscoliq(ig) = absliq1(index,ib) + fint * (absliq1(index+1,ib) - (absliq1(index,ib)))
```

`index` provides a particle size dependent index for the look up table and `ib` provides a waveband dependent index.
The look up table has values for particle sizes in the range [2.5, 59.5] microns in 1 micron intervals (58 values) for each of the 16 different wavebands.
`index` provides a mapping for particle sizes in the range [2.5, 60] microns onto integers in the range [1, 57].
`fint` is between zero and one (except for particle sizes in the range [59.5, 60] microns, in which case it is between 1 and 1.5), and provides a correction term to the absorption coefficient.

For the shortwave, we have:

```fortran
extcoliq(ig) = extliq1(index,ib) + fint * (extliq1(index+1,ib) - extliq1(index,ib))
ssacoliq(ig) = ssaliq1(index,ib) + fint * (ssaliq1(index+1,ib) - ssaliq1(index,ib))
gliq(ig) = asyliq1(index,ib) + fint * (asyliq1(index+1,ib) - asyliq1(index,ib))
forwliq(ig) = gliq(ig)*gliq(ig)
```

The equations for the first three of these shortwave parameters (`extcoliq`, `ssacoliq` and `gliq`, the extinction coefficient, single scattering albedo, and asymmetry parameter respectively) are the same as that for the longwave absorption coefficient, but use different look up tables.
Each of these look up tables has values for particle sizes in the range [2.5, 59.5] microns in 1 micron intervals and for each of the 14 shortwave wavebands.
The forth parameter, `forwliq` is the square of `gliq`.

## Cloud overlap method

This is the RRTMG `icld`.

```python
rrtmg_cloud_overlap_method_dict = {
    'clear_only': 0,
    'random': 1,
    'maximum_random': 2,
    'maximum': 3
}
```

If we choose `clear_only`, there are no clouds, regardless of the other input.
There is no difference between the other three options (`random`, `maximum_random` and `maximum`) for the McICA version of RRTMG.
