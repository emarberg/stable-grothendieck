from tableaux import Tableau


class InsertionState:

    def __init__(self, tableau=None, outer=None):
        self.tableau = Tableau() if tableau is None else tableau
        self.outer = outer
        self.inner_dimensions = self._get_inner_dimensions()

    def _get_inner_dimensions(self):
        if self.outer:
            i, j = self.outer
            t = self.tableau.remove(i, j)
        else:
            t = self.tableau
        m, n = t.max_row(), t.max_column()

        if self.outer:
            i, j = self.outer
            assert {True, False} == {m + 2 == i, n + 2 == j}
        return m, n

    def __eq__(self, other):
        assert type(other) == type(self)
        return self.tableau == other.tableau and self.outer == other.outer

    def __hash__(self):
        return hash((self.tableau, self.outer))

    def __repr__(self):
        return str(self.tableau)

    def __contains__(self, box):
        return box in self.tableau

    def __iter__(self):
        return self.tableau.__iter__()

    @property
    def boxes(self):
        return self.tableau.boxes

    def is_initial(self):
        return self.outer is not None and self.outer[0] == 1

    def is_terminal(self):
        return self.outer is None

    def get(self, i, j, default=None):
        return self.tableau.get(i, j, default)

    def add(self, a):
        assert self.outer is None
        n = self.tableau.max_column() + 2
        tableau = self.tableau.add(1, n, a)
        outer = (1, n)
        return self.__class__(tableau, outer)

    def has_outer_box_in_last_row(self):
        return self.outer and self.outer[0] == self.inner_dimensions[0] + 2

    def has_outer_box_in_last_column(self):
        return self.outer and self.outer[1] == self.inner_dimensions[1] + 2

    def pop(self):
        """
        Returns pair (t, v) where t is the tableau formed by removing the outer box
        and v is value in the outer box. Raises exception if state is terminal.
        """
        assert self.outer is not None
        i, j = self.outer
        return self.tableau.remove(i, j), self.get(i, j)

    def next(self):
        raise NotImplementedError

    def path(self):
        state = self
        while True:
            next_state, bumped = state.next()
            if bumped is None:
                return
            yield next_state, bumped
            state = next_state

    def insert(self, a):
        assert self.is_terminal()
        state = self.add(a)
        path = list(state.path())
        bumping_path = [box for _, box in path]
        final_state = path[-1][0]
        return final_state, bumping_path


class HState(InsertionState):
    pass


class FState(InsertionState):

    @classmethod
    def inverse_insertion(cls, p, q):
        ans = tuple()
        while p:
            q, record = q.pop()
            s, _ = FState(p).previous(record)
            while not s.is_initial():
                s, _ = s.previous()
            p, v = s.pop()
            ans = (v,) + ans
        return ans

    @classmethod
    def insertion_tableaux(cls, *args, **kwargs):
        if len(args) == 1 and type(args[0]) != int:
            args = args[0]
        p = kwargs.get('p', Tableau())
        q = kwargs.get('q', Tableau())
        offset = abs(q.last())
        state = FState(p)
        for index, a in enumerate(args):
            state, path = state.insert(a)
            i, j = path[-1]

            r = index + 1 + offset
            if any(b[0] == b[1] for b in path[:-1]) or ((i, j) not in state and i == j):
                r = -r

            if (i, j) not in state and r > 0:
                a = max([a for (a, b) in state.boxes if b == j - 1])
                i, j = a, j - 1
            elif (i, j) not in state:
                b = max([b for (a, b) in state.boxes if i - 1 == a])
                i, j = i - 1, b

            q = q.add(i, j, r)
        p = state.tableau
        return p, q

    def _row_transition(self):
        i, n = self.outer
        a = self.get(i, n)

        # no row/diagonal transition unless row with outer box excluded is increasing
        row = sorted([self.get(i, x) for x in range(i, n) if (i, x) in self])
        assert row == [self.get(i, x) for x in range(i, i + len(row))]

        # row transition R1
        if all(x <= a for x in row):
            j = i + len(row)
            col = [self.get(r, s) for (r, s) in self.boxes if s == j]
            if all(x < a for x in row) and all(x < a for x in col):
                tableau = self.tableau.add(i, j, a).remove(i, n)
            else:
                tableau = self.tableau.remove(i, n)
            return FState(tableau), (i, j)

        x = i
        while self.get(i, x) <= a:
            x += 1
        b = self.get(i, x)

        if i < x:
            # row transition R2
            if self.get(i, x - 1) < a and (i == 1 or self.get(i - 1, x) < a):
                tableau = self.tableau.add(i + 1, n, b).replace(i, x, a).remove(i, n)
                return FState(tableau, (i + 1, n)), (i, x)
            # diagonal transition D1
            elif x == i + 1 and (i + 1, i + 1) not in self:
                m = self.tableau.max_row()
                tableau = self.tableau.add(m + 2, i + 1, b).remove(i, n)
                return FState(tableau, (m + 2, i + 1)), (i, x)
            # row transition R3
            else:
                tableau = self.tableau.add(i + 1, n, b).remove(i, n)
                return FState(tableau, (i + 1, n)), (i, x)

        if i == x:
            m = self.tableau.max_row()
            # diagonal transition D2
            if a % 2 == b % 2:
                if (i == 1 or self.get(i - 1, i) < a):
                    tableau = self.tableau.replace(i, x, a).add(m + 2, i + 1, b).remove(i, n)
                else:
                    tableau = self.tableau.add(m + 2, i + 1, b).remove(i, n)
            else:
                # raise NotImplementedError
                tableau = self.tableau.add(m + 2, i + 1, b + 1).remove(i, n)
            return FState(tableau, (m + 2, i + 1)), (i, x)

    def _column_transition(self):
        m, j = self.outer
        a = self.get(m, j)

        # no column transition unless column with outer box excluded is increasing
        col = sorted([self.get(x, j) for x in range(1, m) if (x, j) in self])
        assert col == [self.get(x, j) for x in range(1, 1 + len(col))]

        # column transition C1
        if all(x <= a for x in col):
            i = len(col) + 1
            row = [self.get(i, s) for s in range(1, j) if (i, s) in self]
            if all(x < a for x in row) and all(x < a for x in col):
                tableau = self.tableau.add(i, j, a).remove(m, j)
            else:
                tableau = self.tableau.remove(m, j)
            return FState(tableau), (i, j)

        x = 1
        while self.get(x, j) <= a:
            x += 1
        b = self.get(x, j)

        # column transition C2
        if (x == j or self.get(x, j - 1) < a) and (x == 1 or self.get(x - 1, j) < a):
            tableau = self.tableau.add(m, j + 1, b).replace(x, j, a).remove(m, j)
        # column transition C3
        else:
            tableau = self.tableau.add(m, j + 1, b).remove(m, j)
        return FState(tableau, (m, j + 1)), (x, j)

    def next(self):
        if self.has_outer_box_in_last_column():
            return self._row_transition()
        elif self.has_outer_box_in_last_row():
            return self._column_transition()
        else:
            return self, None

    def previous(self, record=None):
        if self.is_initial():
            return self, None
        elif self.has_outer_box_in_last_column():
            return self._inverse_row_transition()
        elif self.has_outer_box_in_last_row():
            return self._inverse_column_transition()
        else:
            return self._inverse_terminal_transition(record)

    def _inverse_terminal_transition(self, record):
        i, j, v = record
        values = {v} if type(v) == int else set(v)
        n = max(abs(x) for x in values)
        assert not (n in values and -n in values)

        if len(values) == 1 and n in values:
            a, b = i, j
            c = self.get(i, j)
            t = self.tableau.remove(i, j)
            j = t.max_column() + 2

        elif len(values) > 1 and n in values:
            t = self.tableau
            while (i - 1, j) in t and (i - 1, j + 1) not in t:
                i = i - 1
            a, b = i, j + 1
            c = max(self.get(i, j), self.get(a - 1, b, 0))
            j = self.tableau.max_column() + 2

        elif len(values) == 1 and -n in values:
            a, b = i, j
            c = self.get(i, j)
            t = self.tableau.remove(i, j)
            i = t.max_row() + 2

        elif len(values) > 1 and -n in values:
            t = self.tableau
            while (i, j - 1) in t and (i + 1, j - 1) not in t and i + 1 <= j - 1:
                j = j - 1
            a, b = i + 1, j
            c = max(self.get(i, j), self.get(a, b - 1, 0))
            i = t.max_row() + 2

        t = t.add(i, j, c)
        return FState(t, (i, j)), (a, b)

    def _inverse_column_transition(self):
        m, j = self.outer
        b = self.get(m, j)

        x = 1
        while (x + 1, j - 1) in self and self.get(x + 1, j - 1) <= b:
            x += 1
        a = self.get(x, j - 1)

        # diagonal transition D1
        if x == j - 1 and b == self.get(j - 1, j) and b % 2 == 0:
            t = self.tableau.remove(m, j)
            i, n = x, t.max_column() + 2
            c = max(t.get(j - 1, j - 1), t.get(j - 2, j, 0))
            t = t.add(i, n, c)
            return FState(t, (i, n)), (x, x)

        # diagonal transition D4
        if x == j - 1 and a < b and b % 2 != 0:
            t = self.tableau.remove(m, j)
            i, n = x, t.max_column() + 2
            t = t.add(i, n, a - 1)
            return FState(t, (i, n)), (x, x)

        # diagonal transition D2
        if x == j - 1 and a == b and self.get(j - 2, j - 1) % 2 == 0:
            t = self.tableau.remove(m, j)
            i, n = x, t.max_column() + 2
            c = t.get(j - 2, j - 1)
            t = t.add(i, n, c)
            return FState(t, (i, n)), (x, x)

        # diagonal transition D3
        if x == j - 1 and a < b and b % 2 == 0:
            t = self.tableau.replace(x, x, b).remove(m, j)
            i, n = x, t.max_column() + 2
            t = t.add(i, n, a)
            return FState(t, (i, n)), (x, x)

        # column transition C3
        if a == b:
            c = max(self.get(x - 1, j - 1), self.get(x, j - 2, 0))
            t = self.tableau.remove(m, j).add(m, j - 1, c)
            return FState(t, (m, j - 1)), (x, j - 1)

        # column transition C4
        else:
            t = self.tableau.remove(m, j).replace(x, j - 1, b).add(m, j - 1, a)
            return FState(t, (m, j - 1)), (x, j - 1)

    def _inverse_row_transition(self):
        i, n = self.outer
        b = self.get(i, n)

        x = i - 1
        while (i - 1, x + 1) in self and self.get(i - 1, x + 1) <= b:
            x += 1
        assert x > i - 1
        a = self.get(i - 1, x)

        # row transition R2
        if a == b:
            c = max(self.get(i - 1, x - 1), self.get(i - 2, x, 0))
            t = self.tableau.remove(i, n).add(i - 1, n, c)
            return FState(t, (i - 1, n)), (i - 1, x)

        # row transition R3
        else:
            t = self.tableau.remove(i, n).replace(i - 1, x, b).add(i - 1, n, a)
            return FState(t, (i - 1, n)), (i - 1, x)
