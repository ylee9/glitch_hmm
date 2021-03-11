function varargout = find_sim_cols(x, varargin)
% FIND_SIM_COLS Computation of indices of columns that are similar in two matrices.
%    [INDX,INDY] = FIND_SIM_COLS(X, Y, THRES) Computes and returns
%    the indices of those columns that are similar, in terms of Euclidean
%    distances, in the two given matrices, X and Y.  Both matrices may
%    have different numbers of columns, but they must have the same number
%    of rows.  The similarity threshold value (of type double) should be
%    specified by the last input parameter, THRES.
%    For instance,
%       x = [1 3 5; 2 4 6];  y = [3 1 4 8 1 5; 1 2.1 7 9 2.11 6.11];
%       [indx,indy] = find_sim_cols(x, y, 0.5);
%    gives
%       indx = [1 1 3]'
%       indy = [2 5 6]'
%    showing that column 1 of x has Euclidean distances smaller than 0.5
%    from colums 2 and 5 of y; column 3 of x has a Euclidean distance
%    smaller than 0.5 from column 6 of y.
%
%    INDX = FIND_SIM_COLS(X, THRES) Computes and returns the indices of
%    those columns in matrix X that are similar, i.e. the Euclidean
%    distances between each other are less than the specified threshold
%    value, THRES (of type double).  The function has been so written
%    that one can later safely delete all the columns marked in INDX.
%    For instance,
%       x = [3 1 4 5 1 5 1; 1 2 7 9 2.11 6 2.1];
%       indx = find_sim_cols(x, 0.5);
%    gives
%       indx = [5 7 7]';
%    Thus,
%       x(:,indx) = [];
%    gives
%       x = [3 1 4 5 5; 1 2 7 9 6].
%
% Sofia Suvorova
% Du Huynh, Oct 2002.

d=[];  siz = size(x,2);
if nargin == 3 & nargout == 2
   y = varargin{1};
   thres = varargin{2};
   for n=1:siz
      d(n,:)=sqrt(sum((y-x(:,n)*ones(1,size(y,2))).^2,1));
   end
   [varargout{1},varargout{2}]=find(d < thres);
elseif nargin == 2 & nargout == 1
   thres = varargin{1};
   for n=1:siz
      d(n,:)=sqrt(sum((x-x(:,n)*ones(1,size(x,2))).^2));
   end
   d = tril(d) + triu(NaN*ones(siz));
   [varargout{1},tmp]=find(d < thres);
else
   error('find_sim_cols: invalid number of arguments');
end
