Parent = [ 2 2 12 5 5 11 8 8 10 10 11 12 -1 ];
Parent = Parent+1;
G = digraph(Parent(Parent~=0),find(Parent));
h = plot(G,'Layout','Layered', "ShowArrows", false);
h.Marker = 'o';
h.MarkerSize = 7;
