clear all
% Calls law_of_cosines.m function, finds sensitivity of measurement error
% in triangle lengths on the accuracy of the triangle angles
% theta(1) is the angle opposite of L(1) etc.
% unknown_guesses: three guesses for the unknowns
% theta: three angles in triangle (radians).  Set any unknown angles to 0
% L: three lengths of the triangle.  Set any unknown lengths to 0
% Any combination of three unknowns are allowed...
%
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

% guesses
unknown=[0.3 2 2];
% theta, radians
theta=[0 0 0];
% L, lengths of triangle sides
L=[57 57 1];

% minmax, sensitivity to +/- this error in length will be found
minmax=0.03125/2;

err_len=-minmax:minmax:minmax;
cnt=0;
for ii=1:length(err_len)
    for jj=1:length(err_len)
        for kk=1:length(err_len)
            cnt=cnt+1;
            A(1)=L(1)+err_len(ii);
            A(2)=L(2)+err_len(jj);
            %A(3)=L(3)+err_len(kk);
            A(3)=L(3);
            unknown=fsolve(@(x) law_of_cosines(x,theta,A),unknown);
            unknown
            tht(cnt)=unknown(2);
            leg(cnt)=A(1);
        end
    end
end

plot(leg,tht*180/pi,'o')
xlabel('Length of side 1')
ylabel('Angle 1, degrees')

unknown