public struct STuple<T1, T2> {
    public readonly T1 Item1;
    public readonly T2 Item2;

    public STuple(T1 a, T2 b) {
        Item1 = a;
        Item2 = b;
    }
}

public struct STuple<T1, T2, T3> {
    public readonly T1 Item1;
    public readonly T2 Item2;
    public readonly T3 Item3;

    public STuple(T1 a, T2 b, T3 c) {
        Item1 = a;
        Item2 = b;
        Item3 = c;
    }
}

public struct STuple<T1, T2, T3, T4> {
    public readonly T1 Item1;
    public readonly T2 Item2;
    public readonly T3 Item3;
    public readonly T4 Item4;

    public STuple(T1 a, T2 b, T3 c, T4 d) {
        Item1 = a;
        Item2 = b;
        Item3 = c;
        Item4 = d;
    }
}

public struct STuple<T1, T2, T3, T4, T5> {
    public readonly T1 Item1;
    public readonly T2 Item2;
    public readonly T3 Item3;
    public readonly T4 Item4;
    public readonly T5 Item5;

    public STuple(T1 a, T2 b, T3 c, T4 d, T5 e) {
        Item1 = a;
        Item2 = b;
        Item3 = c;
        Item4 = d;
        Item5 = e;
    }
}

public struct STuple<T1, T2, T3, T4, T5, T6> {
    public readonly T1 Item1;
    public readonly T2 Item2;
    public readonly T3 Item3;
    public readonly T4 Item4;
    public readonly T5 Item5;
    public readonly T6 Item6;

    public STuple(T1 a, T2 b, T3 c, T4 d, T5 e, T6 f) {
        Item1 = a;
        Item2 = b;
        Item3 = c;
        Item4 = d;
        Item5 = e;
        Item6 = f;
    }
}