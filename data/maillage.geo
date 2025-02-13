l_z_axis = 0.1;
l_bob    = 0.1;
l_inf    = 1;

// parametres de la plaque
e = 1;  // epaisseur de la plaque
L = 20; // longueur de la plaque

// parametres de la bobine
entrefer = 0.2; // hauteur de la bobine au-dessus de la plaque
r_int    = 3; // rayon interieur
r_ext    = 5; // rayon exterieur
h        = 1; // hauteur de la bobine
z_bob    = e/2 + entrefer;

////////////
// plaque //
////////////
p_plaque_1 = newp; Point(p_plaque_1) = {0,-e/2,0,l_z_axis};
p_plaque_2 = newp; Point(p_plaque_2) = {L,-e/2,0,l_inf/4};
p_plaque_3 = newp; Point(p_plaque_3) = {L,e/2,0,l_inf/4};
p_plaque_4 = newp; Point(p_plaque_4) = {0,e/2,0,l_z_axis};

l_plaque_1 = newl; Line(l_plaque_1) = {p_plaque_1,p_plaque_2};
l_plaque_2 = newl; Line(l_plaque_2) = {p_plaque_2,p_plaque_3};
l_plaque_3 = newl; Line(l_plaque_3) = {p_plaque_3,p_plaque_4};
l_plaque_4 = newl; Line(l_plaque_4) = {p_plaque_4,p_plaque_1};

ll_plaque = newll; Line Loop(ll_plaque) = {l_plaque_1,l_plaque_2,l_plaque_3,l_plaque_4};
s_plaque  = news; Surface(s_plaque) = {ll_plaque};

// bobine haute
p_bob_top_1 = newp; Point(p_bob_top_1) = {r_int,z_bob,0,l_bob};
p_bob_top_2 = newp; Point(p_bob_top_2) = {r_ext,z_bob,0,l_bob};
p_bob_top_3 = newp; Point(p_bob_top_3) = {r_ext,z_bob+h,0,l_bob};
p_bob_top_4 = newp; Point(p_bob_top_4) = {r_int,z_bob+h,0,l_bob};

l_bob_top_1 = newll; Line(l_bob_top_1) = {p_bob_top_1,p_bob_top_2};
l_bob_top_2 = newll; Line(l_bob_top_2) = {p_bob_top_2,p_bob_top_3};
l_bob_top_3 = newll; Line(l_bob_top_3) = {p_bob_top_3,p_bob_top_4};
l_bob_top_4 = newll; Line(l_bob_top_4) = {p_bob_top_4,p_bob_top_1};

ll_bob_top = newll; Line Loop(ll_bob_top) = {l_bob_top_1,l_bob_top_2,l_bob_top_3,l_bob_top_4};
s_bob_top  = news; Surface(s_bob_top) = {ll_bob_top};


// bobine basse
p_bob_low_1 = newp; Point(p_bob_low_1) = {r_int,-z_bob,0,l_bob};
p_bob_low_2 = newp; Point(p_bob_low_2) = {r_int,-z_bob-h,0,l_bob};
p_bob_low_3 = newp; Point(p_bob_low_3) = {r_ext,-z_bob-h,0,l_bob};
p_bob_low_4 = newp; Point(p_bob_low_4) = {r_ext,-z_bob,0,l_bob};

l_bob_low_1 = newll; Line(l_bob_low_1) = {p_bob_low_1,p_bob_low_2};
l_bob_low_2 = newll; Line(l_bob_low_2) = {p_bob_low_2,p_bob_low_3};
l_bob_low_3 = newll; Line(l_bob_low_3) = {p_bob_low_3,p_bob_low_4};
l_bob_low_4 = newll; Line(l_bob_low_4) = {p_bob_low_4,p_bob_low_1};

ll_bob_low = newll; Line Loop(ll_bob_low) = {l_bob_low_1,l_bob_low_2,l_bob_low_3,l_bob_low_4};
s_bob_low  = news; Surface(s_bob_low) = {ll_bob_low};

// air haut
p_air_top_1 = newp; Point(p_air_top_1) = {L,L,0,l_inf};
p_air_top_2 = newp; Point(p_air_top_2) = {0,L,0,l_inf};

l_air_top_1 = newl; Line(l_air_top_1) = {p_plaque_3,p_air_top_1};
l_air_top_2 = newl; Line(l_air_top_2) = {p_air_top_1,p_air_top_2};
l_air_top_3 = newl; Line(l_air_top_3) = {p_air_top_2,p_plaque_4};

ll_air_top = newll; Line Loop(ll_air_top) = {-l_plaque_3,l_air_top_1,l_air_top_2,l_air_top_3};
s_air_top = news; Surface(s_air_top) = {ll_air_top,-ll_bob_top};

// air bas
p_air_low_1 = newp; Point(p_air_low_1) = {0,-L,0,l_inf};
p_air_low_2 = newp; Point(p_air_low_2) = {L,-L,0,l_inf};

l_air_low_1 = newl; Line(l_air_low_1) = {p_plaque_1,p_air_low_1};
l_air_low_2 = newl; Line(l_air_low_2) = {p_air_low_1,p_air_low_2};
l_air_low_3 = newl; Line(l_air_low_3) = {p_air_low_2,p_plaque_2};

ll_air_low = newll; Line Loop(ll_air_low) = {l_air_low_1,l_air_low_2,l_air_low_3,-l_plaque_1};
s_air_low  = news; Surface(s_air_low) = {ll_air_low,-ll_bob_low};

///////////////////////
// Physical Surfaces //
///////////////////////
Physical Surface(1) = {s_plaque};
Physical Surface(2) = {s_air_low,s_air_top};
Physical Surface(3) = {s_bob_low};
Physical Surface(4) = {s_bob_top};

Physical Line(1) = {l_plaque_4,l_air_low_1,l_air_low_2,l_air_low_3,l_plaque_2,l_air_top_1,l_air_top_2,l_air_top_3};