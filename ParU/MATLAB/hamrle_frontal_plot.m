Parent = [
];
Parent = Parent+1;
G = digraph(Parent(Parent~=0),find(Parent));
h = plot(G,'Layout','Layered', "ShowArrows", false);