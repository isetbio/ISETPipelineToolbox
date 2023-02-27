function mouseMove (object, eventdata)
% Filler fxn grabbed from MATLAB support team to display the mouse
% coordinates over a figure. Too be expanded upon. 

C = get (gca, 'CurrentPoint');
title(gca, ['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);