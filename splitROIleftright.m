function h = splitROIleftright(h)
h.dat.cl.iscell = double(h.dat.cl.Mrs <h.dat.res.Mrs_thresh & ...
    h.dat.cl.npix <h.dat.cl.npix_high & h.dat.cl.npix >h.dat.cl.npix_low);

% overwrite manual selections
h.dat.cl.iscell(h.dat.cl.manual>1e-3) = 1;
h.dat.cl.iscell(h.dat.cl.manual<-1e-3) = 0;

h.dat.cl.k1 = reshape(h.dat.cl.iscell(h.dat.res.iclust), h.dat.cl.Ly, h.dat.cl.Lx);

