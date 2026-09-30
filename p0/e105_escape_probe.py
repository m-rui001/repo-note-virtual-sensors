r = "\bibitem{a}\tfrac"
s = "\b" + "ibitem{a}" + "\t" + "frac"
print("non-raw literal  bytes:", [b for b in r.encode()][:4], "len:", len(r))
print("chr-built equals:", r == s)
print("shell here-doc passthrough check:")
