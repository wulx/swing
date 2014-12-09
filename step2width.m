function projWidths = step2width(steps, stepAngleDeg, leafWidth, initAngleDeg)
%STEP2WIDTH convert number of steps of the specified stepper motor to
%projected widths

% copyright (c) wulx, gurdy.woo@gmail.com
% last modified by wulx, 2013/10/31

if any(~isfinite(steps))
    error('number of steps should be finite.')
end

if ~isfinite(stepAngleDeg) || stepAngleDeg<0
    error('step angle should be finite non-negative')
end

if nargin<4, initAngleDeg = 0; end

rotAngles = deg2rad(initAngleDeg + stepAngleDeg*steps);

% the scanning leaf parallels to ion beam when rotation angle is 0,
% and perpendiculars to when the angle is  pi/2.
% projWidths will increase as rotAngles increases between 0 and pi/2.
projWidths = abs(sin(rotAngles) * leafWidth);
