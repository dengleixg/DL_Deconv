function ai = funCTAI(p,V50,CNI)

D = p.Dd * (p.b/p.a) * (p.d/p.c);
ai = CNI / (D * (1/(2*V50))^2);