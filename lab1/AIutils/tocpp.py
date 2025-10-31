a = open("multa.txt", "r").readlines()
b = open("multb.txt", "r").readlines()

for i, (A, B) in enumerate(zip(a, b)):
    A = ["+" + x if len(x) == 2 else x for x in A.strip().split(" ")]
    A = [f"{x[0]} A[{int(x[1]) - 1}][{int(x[2]) - 1}]" for x in A]

    B = ["+" + x if len(x) == 2 else x for x in B.strip().split(" ")]
    B = [f"{x[0]} B[{int(x[1]) - 1}][{int(x[2]) - 1}]" for x in B]

    print(f"Matrix h{i + 1} = ({" ".join(A)}) * ({" ".join(B)});")


c = open("multc.txt", "r").readlines()

for i, C in enumerate(c):
    C = ["+" + x if "-" not in x else x for x in C.strip().split(" ")]
    C = [f"{x[0]} h{x[1:]}" for x in C]

    print(f"Matrix c{i % 4 + 1}{i // 4 + 1} = {" ".join(C)};")