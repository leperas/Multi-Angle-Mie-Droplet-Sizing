function F=law_of_cosines(unknown_guesses,theta,L)
% Law of cosines, returns all unknown angles and sides in triangle given 
% the known sides, where theta(1) is the angle opposite of L(1) etc.
% unknown_guesses: three guesses for the unknowns
% theta: three angles in triangle (radians).  Set any unknown angles to 0
% L: three lengths of the triangle.  Set any unknown lengths to 0
% Any combination of three unknowns are allowed...
%
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

unk=0;
for ii=1:3
    if theta(ii)==0
        unk=unk+1;
        theta(ii)=unknown_guesses(unk);
    end
    if L(ii)==0
        unk=unk+1;
        L(ii)=unknown_guesses(unk);
    end
end

F(1)=L(2)^2+L(3)^2-L(1)^2-2*L(2)*L(3)*cos(theta(1));
F(2)=L(3)^2+L(1)^2-L(2)^2-2*L(3)*L(1)*cos(theta(2));
F(3)=L(1)^2+L(2)^2-L(3)^2-2*L(1)*L(2)*cos(theta(3));

end


