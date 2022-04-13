from partitions import Partition
from tableaux import Tableau, FRENCH
from cached import cached_value


VALUED_SET_CACHE = {}
TRANSITION_CACHE = {}
BACKWARD_TRANSITION_CACHE = {}
FORWARD_TRANSITION_CACHE = {}
MIDDLE_TRANSITION_CACHE = {}


class ValuedSetTableau:

    def is_altered(self, index):
        vst = self
        x = vst.hinge(index)
        if x:
            if vst.tableau.get(x, x) == index and vst.tableau.get(x + 1, x + 1) == -index - 1:
                return True
        return False

    def is_separable(self, index):
        h = self.hinge(index)
        if h and self.grouping.get(h, h) == 0 and self.grouping.get(h + 1, h + 1) == 0:
            return True
        return False

    def is_highest_weight(self):
        return all(self.grouping.get(i, j) == 0 for (i, j) in self.tableau.boxes)

    def hinge(self, index):
        assert index > 0
        args = [index, -index, index + 1, -index - 1]
        tab = self.tableau
        x = [x for (x, y) in tab.boxes if x == y and tab.get(x, y) in args and tab.get(x, y + 1) in args and tab.get(x + 1, y + 1) in args]
        assert len(x) <= 1
        return x[0] if x else None

    def is_group_end(self, i, j):
        return (i, j) in self.tableau.boxes and not self.grouping.get(i, j)

    def is_group_start(self, i, j):
        if (i, j) not in self.tableau.boxes:
            return False
        if self.tableau.get(i, j) > 0:
            return (i, j - 1) not in self.tableau.boxes or not (self.tableau.get(i, j - 1) > 0 and self.grouping.get(i, j - 1))
        if self.tableau.get(i, j) < 0:
            return (i + 1, j) not in self.tableau.boxes or not (self.tableau.get(i + 1, j) < 0 and self.grouping.get(i + 1, j))

    def is_singleton(self, i, j):
        return self.is_group_end(i, j) and self.is_group_start(i, j)

    def diagonal_primes(self):
        return len([i for (i, j) in self.tableau.boxes if i == j and self.tableau.get(i, j) < 0])

    def unprime(self):
        t = self.tableau.boxes.copy()
        for x, y in self.singletons():
            t[x, y] = abs(self.tableau.get(x, y))
        return ValuedSetTableau(t, self.grouping)

    def unprime_diagonal(self):
        t = self.tableau.boxes.copy()
        for x in self.diagonal_singletons():
            t[x, x] = abs(self.tableau.get(x, x))
        return ValuedSetTableau(t, self.grouping)

    def singletons(self, *args):
        args = [abs(a) for a in args]
        ans = [(x, y) for (x, y) in self.tableau.boxes if self.is_group_end(x, y) and self.is_group_start(x, y)]
        return [(x, y) for (x, y) in ans if len(args) == 0 or abs(self.tableau.get(x, y)) in args]

    def diagonal_singletons(self, *args):
        return sorted([x for (x, y) in self.singletons(*args) if x == y])

    def primed_diagonal_cells(self, *args):
        args = [-abs(a) for a in args]
        return sorted([x for (x, y) in self.tableau.boxes if x == y and self.tableau.get(x, y) in args])

    def primed_groups(self, *args):
        if len(args) == 0:
            return len([(i, j) for (i, j) in self.tableau.boxes if self.tableau.get(i, j) < 0 and not self.grouping.get(i, j)])
        else:
            args = [-abs(a) for a in args]
            return len([(i, j) for (i, j) in self.tableau.boxes if self.tableau.get(i, j) in args and not self.grouping.get(i, j)])

    def weight(self, n):
        ans = n * [0]
        for i, j in self.grouping.boxes:
            if not self.grouping.get(i, j):
                w = abs(self.tableau.get(i, j))
                ans[w - 1] += 1
        return tuple(ans)

    def __hash__(self):
        return hash((self.tableau, self.grouping))

    def is_semistandard(self, order=None):
        if order is None:
            ofn = (lambda x: -2 - 2 * x if x < 0 else 2 * x - 1)
        if type(order) in [list, tuple]:
            order = {a: i for (i, a) in enumerate(order)}
        if type(order) == dict:
            ofn = (lambda x: order[x])
        rowint = {(i, j) for (i, j) in self.tableau.boxes if (i + 1, j) in self.tableau.boxes}
        colint = {(i, j) for (i, j) in self.tableau.boxes if (i, j + 1) in self.tableau.boxes}
        for (i, j) in rowint:
            p = ofn(self.tableau.get(i, j))
            q = ofn(self.tableau.get(i + 1, j))
            if p > q or (p == q and p % 2 != 0):
                return False
        for (i, j) in colint:
            p = ofn(self.tableau.get(i, j))
            q = ofn(self.tableau.get(i, j + 1))
            if p > q or (p == q and p % 2 == 0):
                return False
        return True

    @classmethod
    def is_valid_position(cls, tab, grp, i, j):
        assert (i, j) in tab.boxes
        grp = grp.get(i, j) if type(grp) != bool else grp
        if not grp:
            return True
        if tab.get(i, j) > 0:
            return tab.get(i, j) == tab.get(i, j + 1)
        else:
            return tab.get(i, j) == tab.get(i - 1, j)

    def __init__(self, entry_mapping, group_mapping):
        self.tableau = entry_mapping if type(entry_mapping) == Tableau else Tableau(entry_mapping)
        self.grouping = group_mapping if type(group_mapping) == Tableau else Tableau(group_mapping)
        assert all(self.is_valid_position(self.tableau, self.grouping, i, j) for (i, j) in self.tableau.boxes)
        self._strval = None

    def get_verticals(self, v):
        assert v < 0
        vertical_starts = sorted([(i, j) for (i, j) in self.tableau.boxes if self.tableau.get(i, j) == v and not self.grouping.get(i, j)], key=lambda ij: (-ij[1], ij[0]))
        vertical_ends = []
        for i, j in vertical_starts:
            while self.grouping.get(i + 1, j) and self.tableau.get(i + 1, j) == self.tableau.get(i, j) < 0:
                i += 1
            vertical_ends.append((i, j))
        return vertical_starts, vertical_ends

    def get_horizontals(self, v):
        assert v > 0
        horizontal_ends = sorted([(i, j) for (i, j) in self.tableau.boxes if self.tableau.get(i, j) == v and not self.grouping.get(i, j)], key=lambda ij: (ij[0], -ij[1]))
        horizontal_starts = []
        for i, j in horizontal_ends:
            while self.grouping.get(i, j - 1) and self.tableau.get(i, j - 1) == self.tableau.get(i, j) > 0:
                j -= 1
            horizontal_starts.append((i, j))
        return horizontal_starts, horizontal_ends

    def forward_transition(self, value):
        return self.cached_forward_transition(self, value)

    @classmethod
    def coincide(cls, t, u, i, j):
        return t.tableau.get(i, j) == u.tableau.get(i, j) and t.grouping.get(i, j) == u.grouping.get(i, j)

    @classmethod
    def undo_alteration(cls, vst, value):
        tab, grp = vst.tableau.boxes.copy(), vst.grouping.boxes.copy()

        t = vst
        x = vst.hinge(value)
        altered = vst.is_altered(value)
        if altered:
            assert vst.tableau.get(x, x) == value and vst.tableau.get(x + 1, x + 1) == -value - 1
            if not vst.grouping.get(x, x) and not vst.grouping.get(x + 1, x + 1):
                tab[x, x] = -value
                t = ValuedSetTableau(Tableau(tab), grp)
                assert not t.is_altered(value)
                u = t.forward_transition(value)
                if cls.coincide(t, u, x, x) and cls.coincide(t, u, x + 1, x + 1):
                    tab[x, x] = value
                    tab[x + 1, x + 1] = value + 1
            elif vst.tableau.get(x, x + 1) < 0 < vst.grouping.get(x + 1, x + 1):
                tab[x, x] = -value
            else:
                tab[x + 1, x + 1] = value + 1
            t = ValuedSetTableau(Tableau(tab), grp)

        return t

    @cached_value(FORWARD_TRANSITION_CACHE)
    def cached_forward_transition(cls, vst, value):  # noqa
        t = cls.undo_alteration(vst, value)

        tab, grp = t.tableau.boxes.copy(), t.grouping.boxes.copy()
        vertical_starts, vertical_ends = t.get_verticals(-value - 1)
        horizontal_starts, horizontal_ends = t.get_horizontals(value)

        for i in range(len(vertical_starts)):
            (vx1, vy1), (vx2, vy2) = vertical_starts[i], vertical_ends[i]
            for j in range(len(horizontal_starts)):
                (hx1, hy1), (hx2, hy2) = horizontal_starts[j], horizontal_ends[j]
                if (vx1, vy1) == (hx1 + 1, hy1):
                    if hy1 < hy2:
                        hy1 += 1
                        vx1 -= 1
                    else:
                        hx1 = vx2
                        hx2 = vx2
                        vx1 -= 1
                        vx2 -= 1
                elif (vx2, vy2) == (hx2, hy2 + 1):
                    if vx1 < vx2:
                        hy2 += 1
                        vx2 -= 1
                    else:
                        vy1 = hy1
                        vy2 = hy1
                        hy1 += 1
                        hy2 += 1
                horizontal_starts[j], horizontal_ends[j] = (hx1, hy1), (hx2, hy2)
            vertical_starts[i], vertical_ends[i] = (vx1, vy1), (vx2, vy2)

        for (vx1, vy1), (vx2, vy2) in zip(vertical_starts, vertical_ends):
            for a in range(vx1, vx2 + 1):
                for b in range(vy1, vy2 + 1):
                    tab[a, b] = -value - 1
                    grp[a, b] = 1
            grp[vx1, vy1] = 0

        for (hx1, hy1), (hx2, hy2) in zip(horizontal_starts, horizontal_ends):
            for a in range(hx1, hx2 + 1):
                for b in range(hy1, hy2 + 1):
                    tab[a, b] = value
                    grp[a, b] = 1
            grp[hx2, hy2] = 0

        return ValuedSetTableau(Tableau(tab), Tableau(grp))

    def middle_transition(self, value, altered, dnp):
        return self.cached_middle_transition(self, value, altered, dnp)

    def _get_row_groups(self, value):
        vst = self
        h1 = [(hx1, hy1, hy2, 0) for ((hx1, hy1), (hx2, hy2)) in zip(*vst.get_horizontals(value))]
        h2 = [(hx1, hy1, hy2, 1) for ((hx1, hy1), (hx2, hy2)) in zip(*vst.get_horizontals(value + 1))]
        hh = sorted(h1 + h2, key=lambda t: (-t[0], t[1]))
        hh_by_rows = []
        for tup in hh:
            if not hh_by_rows or hh_by_rows[-1][-1][0] != tup[0]:
                hh_by_rows.append([])
            hh_by_rows[-1].append(tup)

        one_row_groups = []
        two_row_groups = []

        for i in range(len(hh_by_rows)):
            currow = hh_by_rows[i]
            nextrow = hh_by_rows[i + 1] if i + 1 < len(hh_by_rows) else None
            for (x, y1, y2, f) in currow:
                if nextrow and nextrow[0][0] == x - 1 and nextrow[0][1] <= y2:
                    row2, start2, stop2 = x, y1, currow[-1][2]
                    row1, start1, stop1 = nextrow[0][:3]
                    while hh_by_rows[i + 1]:
                        _, z1, z2, _ = hh_by_rows[i + 1][0]
                        if z1 <= stop2:
                            stop1 = z2
                            hh_by_rows[i + 1] = hh_by_rows[i + 1][1:]
                        else:
                            break
                    two_row_groups.append((row1, start1, stop1, row2, start2, stop2))
                    break
                else:
                    if not one_row_groups or one_row_groups[-1][-1][0][0] != x:
                        one_row_groups.append([0, 0, []])
                    one_row_groups[-1][f] += 1
                    one_row_groups[-1][-1].append((x, y1, y2))

        return one_row_groups, two_row_groups

    def _get_column_groups(self, value):
        vst = self
        v1 = [(vy1, vx1, vx2, 0) for ((vx1, vy1), (vx2, vy2)) in zip(*vst.get_verticals(-value))]
        v2 = [(vy1, vx1, vx2, 1) for ((vx1, vy1), (vx2, vy2)) in zip(*vst.get_verticals(-value - 1))]
        vv = sorted(v1 + v2, key=lambda t: (-t[0], t[1]))
        vv_by_cols = []
        for tup in vv:
            if not vv_by_cols or vv_by_cols[-1][-1][0] != tup[0]:
                vv_by_cols.append([])
            vv_by_cols[-1].append(tup)

        one_col_groups = []
        two_col_groups = []

        for i in range(len(vv_by_cols)):
            currcol = vv_by_cols[i]
            nextcol = vv_by_cols[i + 1] if i + 1 < len(vv_by_cols) else None
            for (y, x1, x2, f) in currcol:
                if nextcol and nextcol[0][0] == y - 1 and nextcol[0][1] <= x2:
                    col2, start2, stop2 = y, x1, currcol[-1][2]
                    col1, start1, stop1 = nextcol[0][:3]
                    while vv_by_cols[i + 1]:
                        _, z1, z2, _ = vv_by_cols[i + 1][0]
                        if z1 <= stop2:
                            stop1 = z2
                            vv_by_cols[i + 1] = vv_by_cols[i + 1][1:]
                        else:
                            break
                    two_col_groups.append((col1, start1, stop1, col2, start2, stop2))
                    break
                else:
                    if not one_col_groups or one_col_groups[-1][-1][0][0] != y:
                        one_col_groups.append([0, 0, []])
                    one_col_groups[-1][f] += 1
                    one_col_groups[-1][-1].append((y, x1, x2))

        return one_col_groups, two_col_groups

    @cached_value(MIDDLE_TRANSITION_CACHE)
    def cached_middle_transition(cls, vst, value, altered, dnp):  # noqa
        tab, grp = vst.tableau.boxes.copy(), vst.grouping.boxes.copy()
        case = '*' if altered else None

        one_row_groups, two_row_groups = vst._get_row_groups(value)
        one_col_groups, two_col_groups = vst._get_column_groups(value)

        for p, q, g in one_row_groups:
            assert len(g) == p + q
            for i in range(len(g)):
                row, col1, col2 = g[i]
                for y in range(col1, col2 + 1):
                    tab[row, y] = (value + 1) if i < q else value
                    grp[row, y] = 1
                grp[row, col2] = 0

        for row1, start1, stop1, row2, start2, stop2 in two_row_groups:
            for y in range(start1, stop1 + 1):
                tab[row1, y] = value + 1
            for y in range(start2, stop2 + 1):
                tab[row2, y] = value
            for y in range(max(start1, start2), min(stop1, stop2)):
                grp[row1, y], grp[row2, y] = grp[row2, y], grp[row1, y]
            # special diagonal condition
            if row1 == start1:
                x = row1
                if not altered and vst.is_group_end(x, x):
                    if dnp and vst.grouping.get(x, x + 1) and vst.grouping.get(x + 1, x + 2) is not None:
                        case = 'a1'
                        grp[x + 1, x + 1] = 0
                        grp[x, x] = 1
                    else:
                        case = 'a2'
                        tab[x, x] = -value
                elif altered and vst.is_group_end(x + 1, x + 1) and vst.grouping.get(x + 1, x + 2) is not None:
                    case = 'a3'
                    grp[x, x] = 0
                    grp[x, x + 1] = 1
                elif altered:
                    assert not vst.is_group_end(x, x)
                    case = 'a4'

        for p, q, g in one_col_groups:
            assert len(g) == p + q
            for i in range(len(g)):
                col, row1, row2 = g[i]
                for x in range(row1, row2 + 1):
                    tab[x, col] = (-value - 1) if i < q else -value
                    grp[x, col] = 1
                grp[row1, col] = 0

        for col1, start1, stop1, col2, start2, stop2 in two_col_groups:
            for x in range(start1, stop1 + 1):
                tab[x, col1] = -value - 1
            for x in range(start2, stop2 + 1):
                tab[x, col2] = -value
            for x in range(max(start1, start2) + 1, min(stop1, stop2) + 1):
                grp[x, col1], grp[x, col2] = grp[x, col2], grp[x, col1]
            # special diagonal condition
            if col2 == stop2:
                x = col2
                if not altered and vst.is_group_end(x, x):
                    if dnp and vst.grouping.get(x - 1, x) and vst.grouping.get(x - 2, x - 1) is not None:
                        case = 'b1'
                        grp[x - 1, x - 1] = 0
                        grp[x, x] = 1
                    else:
                        case = 'b2'
                        tab[x, x] = value + 1
                elif altered and vst.is_group_end(x - 1, x - 1) and vst.grouping.get(x - 2, x - 1) is not None:
                    case = 'b3'
                    grp[x, x] = 0
                    grp[x - 1, x] = 1
                elif altered:
                    assert not vst.is_group_end(x, x)
                    case = 'b4'

        ans = ValuedSetTableau(Tableau(tab), Tableau(grp))
        h = vst.hinge(value)
        if case == '*' and ans.tableau.get(h, h) < 0 and ans.tableau.get(h, h + 1) < 0 and ans.tableau.get(h + 1, h + 1) > 0:
             ans = ValuedSetTableau(ans.tableau.set(h + 1, h + 1, ans.tableau.get(h + 1, h + 1) * -1), ans.grouping)
             case = 'b0'
        elif case == '*' and ans.tableau.get(h + 1, h + 1) > 0 and ans.tableau.get(h, h + 1) > 0 and ans.tableau.get(h, h) < 0:
             ans = ValuedSetTableau(ans.tableau.set(h, h, ans.tableau.get(h, h) * -1), ans.grouping)
             case = 'a0'
        assert case != '*'

        return ans, case

    def backward_transition(self, value):
        return self.cached_backward_transition(self, value)

    @cached_value(BACKWARD_TRANSITION_CACHE)
    def cached_backward_transition(cls, vst, value):  # noqa
        tab, grp = vst.tableau.boxes.copy(), vst.grouping.boxes.copy()
        vertical_starts, vertical_ends = vst.get_verticals(-value)
        horizontal_starts, horizontal_ends = vst.get_horizontals(value + 1)

        for i in range(len(vertical_starts) - 1, -1, -1):
            (vx1, vy1), (vx2, vy2) = vertical_starts[i], vertical_ends[i]
            for j in range(len(horizontal_starts) - 1, -1, -1):
                (hx1, hy1), (hx2, hy2) = horizontal_starts[j], horizontal_ends[j]
                if (vx1, vy1) == (hx1, hy1 - 1):
                    if vx1 < vx2:
                        vx1 += 1
                        hy1 -= 1
                    else:
                        vy1 = hy2
                        vy2 = hy2
                        hy1 -= 1
                        hy2 -= 1
                elif (vx2, vy2) == (hx2 - 1, hy2):
                    if hy1 < hy2:
                        vx2 += 1
                        hy2 -= 1
                    else:
                        hx1 = vx1
                        hx2 = vx1
                        vx1 += 1
                        vx2 += 1
                horizontal_starts[j], horizontal_ends[j] = (hx1, hy1), (hx2, hy2)
            vertical_starts[i], vertical_ends[i] = (vx1, vy1), (vx2, vy2)

        for (vx1, vy1), (vx2, vy2) in zip(vertical_starts, vertical_ends):
            for a in range(vx1, vx2 + 1):
                for b in range(vy1, vy2 + 1):
                    tab[a, b] = -value
                    grp[a, b] = 1
            grp[vx1, vy1] = 0

        for (hx1, hy1), (hx2, hy2) in zip(horizontal_starts, horizontal_ends):
            for a in range(hx1, hx2 + 1):
                for b in range(hy1, hy2 + 1):
                    tab[a, b] = value + 1
                    grp[a, b] = 1
            grp[hx2, hy2] = 0

        return ValuedSetTableau(Tableau(tab), Tableau(grp))

    def transition(self, index, dnp):
        return self.cached_transition(self, index, dnp)

    @classmethod
    def p_adjust(cls, ans, index):
        tab = ans.tableau
        grp = ans.grouping
        h = ans.hinge(index)

        if h:
            if ans.is_singleton(h, h):
                tab = tab.set(h, h, abs(tab.get(h, h)))
            if ans.is_singleton(h + 1, h + 1):
                tab = tab.set(h + 1, h + 1, abs(tab.get(h + 1, h + 1)))

        elif ans.primed_diagonal_cells(index, index + 1):
            x = y = list(ans.primed_diagonal_cells(index, index + 1))[0]
            if tab.get(x, y) == -index - 1:
                tab = tab.set(x, x, index + 1)
                while tab.get(x, y) == index + 1:
                    y += 1
                while grp.get(x, y):
                    y += 1
                for z in range(x, y):
                    grp = grp.set(x, z, grp.get(x, z + 1))
                    tab = tab.set(x, z, index + 1)
                grp = grp.set(x, y, 0)
                tab = tab.set(x, y, -index)

            elif tab.get(x, y) == -index:
                while tab.get(x, y) == -index:
                    x -= 1
                assert grp.get(x, y) == 0

                for z in range(x, y):
                    grp = grp.set(z, y, grp.get(z + 1, y))
                    tab = tab.set(z, y, tab.get(z + 1, y))
                grp = grp.set(y, y, 0)
                tab = tab.set(y, y, index)

                tab = tab.set(x, y, -index - 1)
                while grp.get(x + 1, y) != 0:
                    x += 1
                    tab = tab.set(x, y, -index - 1)

        return ValuedSetTableau(tab, grp)

    @classmethod
    def q_adjust(cls, ans, index, h, case):
        if h:
            tab = ans.tableau
            grp = ans.grouping
            if case in ['a1', 'a2', 'a4']:
                if ans.is_singleton(h + 1, h + 1):
                    tab = tab.set(h + 1, h + 1, tab.get(h + 1, h + 1) * -1)
            if case in ['b1', 'b2', 'b4']:
                if ans.is_singleton(h, h):
                    tab = tab.set(h, h, tab.get(h, h) * -1)
            ans = ValuedSetTableau(tab, grp)
        return ans

    @classmethod
    def reorient(cls, ans, index):
        mapping = {index: index + 1, -index: -index - 1, index + 1: index, -index - 1: -index}
        boxes = {}
        for i, j in ans.tableau.boxes:
            v = ans.tableau.get(i, j)
            boxes[i, j] = mapping.get(v, v)
        return ValuedSetTableau(Tableau(boxes), ans.grouping)

    @cached_value(TRANSITION_CACHE)
    def cached_transition(cls, vst, index, dnp):  # noqa
        altered = vst.is_altered(index)
        f = vst.forward_transition(index)
        m, case = f.middle_transition(index, altered, dnp)
        ans = m.backward_transition(index)

        if not dnp:
            ans = cls.p_adjust(ans, index)
        else:
            h = vst.hinge(index)
            ans = cls.q_adjust(ans, index, h, case)
        return cls.reorient(ans, index)


    @classmethod
    def all(cls, max_entry, mu, nu=(), diagonal_nonprimes=True):
        return cls._all(max_entry, mu, nu, diagonal_nonprimes)

    @cached_value(VALUED_SET_CACHE)
    def _all(cls, max_entry, mu, nu, diagonal_nonprimes):  # noqa
        ans = set()
        tabs = Tableau.semistandard_shifted(max_entry, mu, nu, diagonal_nonprimes)
        for tab in tabs:
            grp = {(i, j): int(cls.is_valid_position(tab, True, i, j)) for (i, j) in tab.boxes}
            g = sorted(grp)
            e = sum(grp.values())
            for v in range(2**e):
                next_grp = {}
                for i, j in g:
                    if grp[i, j]:
                        next_grp[i, j] = v % 2
                        v = v // 2
                    else:
                        next_grp[i, j] = 0
                next_grp = Tableau(next_grp)
                ans.add(cls(tab, next_grp))
        return ans

    def compute_string_array(self):
        boxes = {
            k: ','.join([
                str(-x) + '\'' if x < 0 else str(x)
                for x in sorted(v, key=lambda x:(abs(x), x))
            ]) for k, v in self.tableau.boxes.items()}

        allmax = max([2] + [len(str(boxes[b])) for b in boxes])
        column_spacing_offset = 0

        def pad(x, j):
            return str(x) + (column_spacing_offset + allmax - len(str(x))) * ' '

        array = []
        for i in range(1, self.tableau.max_row() + 1):
            row = []
            for j in range(1, self.tableau.max_column() + 1):
                row += [pad(boxes.get((i, j), '.'), j)]
            array += [row]
        return array, allmax

    def in_same_group(self, i, j, x, y):
        if i == x and y == j + 1:
            return bool(self.grouping.get(i, j, default=False)) and self.tableau.get(i, j) > 0
        if i == x and y + 1 == j:
            return self.in_same_group(x, y, i, j)
        if i + 1 == x and y == j:
            return bool(self.grouping.get(i, j, default=False)) and self.tableau.get(i, j) < 0
        if i == x + 1 and y == j:
            return self.in_same_group(x, y, i, j)
        return False

    def __lt__(self, other):
        return self.tableau < other.tableau and self.grouping < other.grouping

    def __eq__(self, other):
        return type(other) == type(self) and self.tableau == other.tableau and self.grouping == other.grouping

    def __repr__(self):
        if self._strval is None:
            array, column_width = self.compute_string_array()
            m, n = len(array), (len(array[0]) if array else 0)
            boxes = self.tableau.boxes
            newarray = [(2 * n + 1) * [" "] for _ in range(2 * m + 1)]

            for i in range(2 * m + 1):
                for j in range(2 * n + 1):
                    a, b = i // 2 + 1, j // 2 + 1
                    if i % 2 != 0 and j % 2 != 0:
                        newarray[i][j] = array[i // 2][j // 2]
                    elif i % 2 == 0 and j % 2 != 0:
                        bordering = (a, b) in boxes or (a - 1, b) in boxes
                        newarray[i][j] = column_width * ('─' if bordering else ' ')
                    elif i % 2 != 0 and j % 2 == 0:
                        bordering = (a, b) in boxes or (a, b - 1) in boxes
                        newarray[i][j] = '│' if bordering else ' '
                    else:
                        if ((a - 1, b) in boxes and (a, b - 1) in boxes) or ((a, b) in boxes and (a - 1, b - 1) in boxes):
                            newarray[i][j] = '┼'
                        elif (a, b) in boxes and (a, b - 1) in boxes:
                            newarray[i][j] = '┴'
                        elif (a - 1, b) in boxes and (a - 1, b - 1) in boxes:
                            newarray[i][j] = '┬'
                        elif (a, b) in boxes and (a - 1, b) in boxes:
                            newarray[i][j] = '├'
                        elif (a, b - 1) in boxes and (a - 1, b - 1) in boxes:
                            newarray[i][j] = '┤'
                        elif (a, b) in boxes:
                            newarray[i][j] = '└'
                        elif (a - 1, b) in boxes:
                            newarray[i][j] = '┌'
                        elif (a, b - 1) in boxes:
                            newarray[i][j] = '┘'
                        elif (a - 1, b - 1) in boxes:
                            newarray[i][j] = '┐'

            array = newarray
            for (a, b) in boxes:
                i, j = 2 * a - 1, 2 * b - 1
                if not self.grouping.get(a, b):
                    continue
                if self.tableau.get(a, b) > 0:
                    array[i][j + 1] = ' '
                    c = array[i + 1][j + 1]
                    array[i + 1][j + 1] = '┴' if c == '┼' else '─'
                    c = array[i - 1][j + 1]
                    array[i - 1][j + 1] = '┬' if c == '┼' else '─'
                else:
                    array[i - 1][j] = ' ' * column_width
                    c = array[i - 1][j + 1]
                    array[i - 1][j + 1] = '├' if c == '┼' else '│'
                    c = array[i - 1][j - 1]
                    array[i - 1][j - 1] = '┤' if c == '┼' else '│'

            array = array if not FRENCH else list(reversed(array))
            self._strval = '\n\n' + '\n'.join([''.join(line) for line in array]) + '\n\n'
        return self._strval
