function [stens,segment] = GenerateStencils(x,ns,n)
   N = length(x);
   stens = zeros(N,n);
   segment = zeros(N,ns);
   for i=1:N
      x0 = x(i); 
      r = sqrt((x(:)-x0).^2 );
      [~,ix] = sort(r);
      stens(i,1:n) = ix(1:n);
      segment(i,1:ns) = ix(1:ns);
      % include the i point
   end % for i
end