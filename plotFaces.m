% function plotFaces(nel,x,elsize,volfrac)
% figure(1)
% clf; caxis([0,1]); colormap hot; set(gcf,'GraphicsSmoothing','off');
% %% v is x with void margins
% v = zeros(elsize+2,'single');
% s1 = elsize+1;
% v(2:s1(1),2:s1(2),2:s1(3)) = reshape(x,elsize);
% %% x0 is cut level to plot only volfrac part
% n = ceil(volfrac*nel);
% s = sort(x,'descend');
% x0 = (s(n) + s(n+1))/2;
% %% Find out number of boundary faces and reserve memory
% s = (v < x0);
% n = nnz(diff(s,1,1)) + nnz(diff(s,1,2)) + nnz(diff(s,1,3));
% vx = zeros(n,4,'int32');     % x coordinates
% vy = zeros(n,4,'int32');     % y coordinates
% vz = zeros(n,4,'int32');     % z coordinates
% col = zeros(n,1,'single');   % colors
% %% function to add faces
% n = 0;
% d = int32([[0 0 0 0]; [0 0 1 1]; [0 1 1 0]; [1 1 1 1];]);
%     function addFace(dj,di,dk);
%         n = n+1;
%         vx(n,:) = i+d(di,:);
%         vy(n,:) = j+d(dj,:);
%         vz(n,:) = k+d(dk,:);
%         col(n) = v(j,i,k);
%     end
% %% Scan neighbors of 'solid' elements to detect boundary faces
% for k = 2:s1(3);  for i = 2:s1(2);  for j = 2:s1(1);
%             if ~s(j,i,k)
%                 if s(j,i,k-1) addFace(2,3,1); end
%                 if s(j,i,k+1) addFace(2,3,4); end
%                 if s(j,i-1,k) addFace(2,1,3); end
%                 if s(j,i+1,k) addFace(2,4,3); end
%                 if s(j-1,i,k) addFace(1,3,2); end
%                 if s(j+1,i,k) addFace(4,3,2); end
%             end
%         end;  end;  end;
% %% Draw
% patch(vx',vz',-vy',col);
% light('Position',[1 0 1]); view([30,30]);
% axis equal; axis tight; axis off; drawnow;
% end
function plotFaces(nel,x,elesize,volfrac,symmdraw)
% symmdraw=1;
clf; caxis([0,1]); colormap hot; set(gcf,'GraphicsSmoothing','off');
	%% v is x with void margins
	v = zeros(elesize+2,'single');
	s1 = elesize+1;
	v(2:s1(1),2:s1(2),2:s1(3)) = reshape(x,elesize);
	%% x0 is cut level to plot only volfrac part
	n = ceil(volfrac*nel);
	s = sort(x,'descend');
  x0 = (s(n) + s(n+1))/2;
	%% Find out number of boundary faces and reserve memory
  s = (v < x0);
	n = nnz(diff(s,1,1)) + nnz(diff(s,1,2)) + nnz(diff(s,1,3));
  vx = zeros(n,4,'int32');     % x coordinates
	vy = zeros(n,4,'int32');     % y coordinates
	vz = zeros(n,4,'int32');     % z coordinates
	col = zeros(n,1,'single');   % colors
	%% function to add faces
  n = 0;
  d = int32([[0 0 0 0]; [0 0 1 1]; [0 1 1 0]; [1 1 1 1];]);
	function addFace(dj,di,dk);
		n = n+1; 
		vx(n,:) = i+d(di,:);
		vy(n,:) = j+d(dj,:);
		vz(n,:) = k+d(dk,:);
		col(n) = v(j,i,k);
  end
  %% Scan neighbors of 'solid' elements to detect boundary faces
  for k = 2:s1(3);  for i = 2:s1(2);  for j = 2:s1(1);
		if ~s(j,i,k)
			if s(j,i,k-1) addFace(2,3,1); end
			if s(j,i,k+1) addFace(2,3,4); end
			if s(j,i-1,k) addFace(2,1,3); end
			if s(j,i+1,k) addFace(2,4,3); end
			if s(j-1,i,k) addFace(1,3,2); end
			if s(j+1,i,k) addFace(4,3,2); end
		end
	end;  end;  end;
	%% Draw
  patch(vx',vz',-vy',col);
	if symmdraw
		patch(vx',2*elesize(3)+4-vz',-vy',col); % symmetric draw for negative y and z
	end;
  light('Position',[1 0 1]); view([30,30]); 
	axis equal; axis tight; axis off; drawnow;
end