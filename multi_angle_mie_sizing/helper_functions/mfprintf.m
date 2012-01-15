function [count]=mfprintf(fids,varargin)
% MFPRINTF Multiple FID version of writing formatted 
%          data to file.
%    COUNT = MFPRINTF(FIDS,FORMAT,A,...) 
%    Behavior is identical with FPRINTF if FIDS is only 
%    a single FID.  See HELP FPRINTF for details.
%
%    If FIDS is a 1-D array of values the same output goes to 
%    each FID location listed in the FIDS array.
%
%    Example:
%    shoe_length=7.1;
%    fids=[1 3 4];
%    [count]=mfprintf(fids,'Length of shoe: %f \n', shoe_length);
%
%    Output into screen and into FID 3 and 4:
%       Length of shoe: 7.100000
%
% Copyright (C) 2012, Stephen D. LePera.  GNU General Public License, v3.
% Terms available in GPL_v3_license.txt, or <http://www.gnu.org/licenses/>

for ii=1:numel(fids)
    [count(ii)]=fprintf(fids(ii),varargin{:});
end
end
