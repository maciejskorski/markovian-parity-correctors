(* Wolfram check of scripts/prove_rate.py: sharp exponent of max|bias(S_n)| over binary chains with next-bit bias <= beta.
   d1=(1+b u)/2, d2=(1+b v)/2, b=beta in [0,1], u,v in [-1,1]; K = P.diag(1,-1).
   Columns (e = ceil(n/2)): bound |bias|<=b^e (ForAll, exact real QE); b^e attained for some 0<b<1 (Exists); same for b^2.
   Expected: n=3,4: b^2 = b^e is sharp; n=5,6: b^2 is NOT attained (False), the sharp exponent is e=3. *)
K = {{(1 + b u)/2, (1 + b u)/2 - 1}, {1 - (1 + b v)/2, -(1 + b v)/2}};
box = 0 <= b <= 1 && -1 <= u <= 1 && -1 <= v <= 1;
attained[bias_, e_] := Reduce[Exists[{b, u, v}, box && 0 < b < 1 && Abs[bias] == b^e], Reals];
Grid[Prepend[Flatten[Table[
    With[{bias = (MatrixPower[K, n] . {1, 1})[[x0 + 1]], e = Ceiling[n/2]},
      {n, x0, Reduce[ForAll[{b, u, v}, box, Abs[bias] <= b^e], Reals], attained[bias, e], attained[bias, 2]}],
    {n, {3, 4, 5, 6}}, {x0, {0, 1}}], 1],
  {"n", "x0", "|bias|<=b^ceil(n/2)", "b^ceil(n/2) attained", "b^2 attained"}]]
