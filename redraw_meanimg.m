
function redraw_meanimg(h)

I = h.dat.mimg(:,:,h.dat.map);
Irange = I(max(1,h.dat.ylim(1)):min(size(I,1),h.dat.ylim(2)), ...
    max(1,h.dat.xlim(1)):min(size(I,1),h.dat.xlim(2)));

axes(h.axes2); imagesc(I, [min(Irange(:)) max(Irange(:))]);
xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off
colormap('gray')
axes(h.axes3); imagesc(I, [min(Irange(:)) max(Irange(:))]);
xlim([h.dat.xlim]); ylim([h.dat.ylim]);
axis off
colormap('gray')