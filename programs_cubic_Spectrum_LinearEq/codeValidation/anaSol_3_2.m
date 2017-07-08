anaSol_3;
% use robin boundary condition to simulate second type boundary condition
problemPars.boundaryType='robin';
problemPars.robinF_ph={{ {@(y)y,@(z)-1} ; 0 ; {@(y)y,@(z)-z} };   % [-y;0;-yz]
                       { {@(y)y,@(z) 1} ; 0 ; {@(y)y,@(z) z} };   % [ y;0; yz]
                       { {@(z)-1,@(x)x} ; 0 ; {@(z)-z,@(x)x} };   % [-x;0;-xz]
                       { {@(z) 1,@(x)x} ; 0 ; {@(z) z,@(x)x} };   % [ x;0; xz]
                       { 0 ; 0 ; {@(x)-x,@(y)y} };                % [ 0;0;-xy]
                       { 0 ; 0 ; {@(x) x,@(y)y} } };              % [ 0;0; xy]
problemPars.robinA_ph=0;