function setWindowPosition(fig,vPix,hPix)
% Make figure display at the middle of the screen
ss = get(0,'screensize'); %The screen size
width = ss(3);
height = ss(4);

set(fig,'Position',[(width/2)-hPix/2, (height/2)-vPix/2, hPix vPix]);

end