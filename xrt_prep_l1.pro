;#NAME: xrt_prep_l1
;
;#PURPOSE: to generate calibrated level 1 data from level 0 XRT data, including
;         despike and removing contaminated spots on CCD, and dejitter.

;#CALLING SEQUENCE: xrt_prep_l1, dir_xrt
;
;#KEYWORDS: dir_xrt: directory containing XRT level0 fits files
;
;#OUTPUT: level 1 fits files of XRT in a new directory '/l1/'
;
;
;#AUTHOR: Zihao Yang, Peking Univerisity, July 12th, 2019, at UCL-MSSL.


pro xrt_prep_l1, dir_xrt

ff=file_search('./*XRT*.fits')
n=n_elements(ff)

read_xrt, ff, index0, data0
xrt_prep, index0, data0, index1, data1, /despike_despot, coalign=2, /normalize,/float,/miss_map
xrt_jitter, index1, off
data_ca=shift_img(data1, off)

for ii=0, n-1 do begin
  write_xrt,index1[ii],data_ca[*,*,ii],outdir=dir_xrt,outfile=strmid(ff[ii],24,21,/rev)+'l1.fits'
end

outdir='./l1/'
file_mkdir, outdir
file_move,'*.l1.fits', outdir

end
