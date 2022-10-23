n=0;
VM=3000.;   % velocity. ft/sec
DEL=5./57.3;    % delta deflection. radians. 
ALT=0.;         % altitude. ft. 
A=1000.;            % speed of sound 
DIAM=1.;        % diameter. feet.
FR=3.;
XL=20.;
CTW=0.;
CRW=6.;
HW=2.;  
CTT=0.;  % tail center 
CRT=2.; % tail length
HT=2.;
XN=4.;
XCG=10.;% ft
XHL=19.5;
WGT=1000.;% weight 1000lb. = 453.592 kg. 
if ALT<=30000.
    RHO=.002378*exp(-ALT/30000.); % slug/ft^33, 0.002378 slugs/ft^3=1.226kg/m^3
else
    RHO=.0034*exp(-ALT/22000.);
end
SWING=.5*HW*(CTW+CRW); % wind area
STAIL=.5*HT*(CTT+CRT);  % tail area 
SREF=3.1416*DIAM*DIAM/4.; % reference area. ft^2. 
XLP=FR*DIAM;
SPLAN=(XL-XLP)*DIAM+1.33*XLP*DIAM/2.; % planform area of c cylendrical body with a parabolic nose (radome). approximation. 
XCPN=2*XLP/3;   % nose length 
AN=.67*XLP*DIAM;    % nose area
AB=(XL-XLP)*DIAM;   % body area 
XCPB=(.67*AN*XLP+AB*(XLP+.5*(XL-XLP)))/(AN+AB); % body center of pressure
XCPW=XLP+XN+.7*CRW-.2*CTW;  % wing center of pressure
XMACH=VM/A; % mach number 
XIYY=WGT*(3*((DIAM/2)^2)+XL*XL)/(12*32.2);  % moment of inertia. lbf*ft/s^2 
TMP1=(XCG-XCPW)/DIAM;
TMP2=(XCG-XHL)/DIAM;
TMP3=(XCG-XCPB)/DIAM;
TMP4=(XCG-XCPN)/DIAM;

B=sqrt(XMACH^2-1);  % normalized speed for supersonic travel. 
Q=.5*RHO*VM*VM;     % 10701 lb/ft^2. dynamic pressure. lb/ft^2
THD=0;  %   q. theta dot
ALF=0;  % angle of attack. alpha
T=0;

H=.0025;
S=0.;

while T<1.99999
    THDOLD=THD;
    ALFOLD=ALF;
    STEP=1;
    FLAG=0;
    while STEP<=1
        if FLAG==1
            STEP=2;
            THD=THD+H*THDD;
            ALF=ALF+H*ALFD;
            T=T+H;
        end
        % normal force coefficient:
        CN=2*ALF+1.5*SPLAN*ALF*ALF/SREF+8*SWING*ALF/(B*SREF)+...
                8*STAIL*(ALF+DEL)/(B*SREF);
        % pitch moment coefficient:
        CM=2*ALF*TMP4+1.5*SPLAN*ALF*ALF*TMP3/SREF+...
                8*SWING*ALF*TMP1/(B*SREF)...
                +8*STAIL*(ALF+DEL)*TMP2/(B*SREF);
        THDD=Q*SREF*DIAM*CM/XIYY; % q_dot = Qs*d*Cm / i_yy
        XNL=32.2*Q*SREF*CN/WGT;
        ALFD=THD-XNL/VM;
        FLAG=1;
    end

    FLAG=0;
    THD=.5*(THDOLD+THD+H*THDD);
    ALF=.5*(ALFOLD+ALF+H*ALFD);
    S=S+H;
    if S>=.0099999
        S=0.;
        n=n+1;
        ArrayT(n)=T;
        ArrayXNLG(n)=XNL/32.2;
        ArrayALFDEG(n)=ALF*57.3;
    end
end

figure
plot(ArrayT,ArrayXNLG),grid
xlabel('Time (Sec)')
ylabel('Missile Acceleration (G)')

figure
plot(ArrayT,ArrayALFDEG),grid
xlabel('Time (Sec)')
ylabel('Angle of Attack (Deg)')


% clc
output=[ArrayT',ArrayXNLG',ArrayALFDEG'];
save datfil.txt output -ascii
disp 'simulation finished'




dlta = (-15 : 15) * pi / 180;
alf = (-20 : 20) * pi / 180;

CN = zeros(length(dlta), length(alf));
CM = zeros(length(dlta), length(alf));

for d = 1 : length(dlta)
    for a = 1 : length(alf)
        CN(d, a) = 2*alf(a) + 1.5*SPLAN*alf(a)*alf(a)/SREF + 8*SWING*alf(a)/(B*SREF) ...
                    + 8*STAIL*(alf(a)+dlta(d))/(B*SREF);
        CM(d, a) = 2*alf(a)*TMP4 + 1.5*SPLAN*alf(a)*alf(a)*TMP3/SREF ...
                    + 8*SWING*alf(a)*TMP1/(B*SREF)...
                    + 8*STAIL*(alf(a)+dlta(d))*TMP2/(B*SREF);
    end
end


figure(1)
hold on 
figure(2)
hold on 

for d = 1 : 3 : length(dlta)
    c = rand(1, 3);
    figure(1)
    plot(alf * 180 / pi, CN(d, :), 'color', c)
    figure(2)
    plot(alf * 180 / pi, CM(d, :), 'color', c)
end

figure(1)
hold off
grid
box
title('Normal Force Coefficient', 'fontname', 'times', 'fontsize', 12)
xlabel('angle of attack (deg)', 'fontname', 'times', 'fontsize', 12)
ylabel('C_N', 'fontname', 'times', 'fontsize', 12)
lg = legend(strsplit(num2str(dlta(1 : 3 : end)*180/pi)), 'location', 'best');
lg_ttl = get(lg, 'title');
set(lg_ttl, 'string', '\delta')

figure(2)
hold off
grid
box
title('Pitch Moment Coefficient', 'fontname', 'times', 'fontsize', 12)
xlabel('angle of attack (deg)', 'fontname', 'times', 'fontsize', 12)
ylabel('C_M', 'fontname', 'times', 'fontsize', 12)
lg = legend(strsplit(num2str(dlta(1 : 3 : end)*180/pi)), 'location', 'best');
lg_ttl = get(lg, 'title');
set(lg_ttl, 'string', '\delta')

z_db = struct('Alpha', alf, 'Delta', dlta, 'Cz', CN, 'Cm', CM ...
                , 'v', VM * 0.3048, 'diam', DIAM * 0.3048 ...
                , 'cg', XCG * 0.3048, 'iyy', XIYY * 0.04214011 ...
                , 's', SREF / (1/0.3048)^2);
save('aerodynamic_coef', 'z_db')
saveas(1, 'Cz.jpg')
saveas(2, 'Cm.jpg')







