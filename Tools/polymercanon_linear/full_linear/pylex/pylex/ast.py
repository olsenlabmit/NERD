"""Abstract syntax tree class."""

from .nfa import NFA, NFAState


class AST:
    """Node in a regular expression abstract syntax tree."""

    def to_nfa(self, accepting_id=1):
        """Convert this AST to an NFA.

        Arguments:
        accepting_id -- The ID of the accepting state. Defaults to 1.

        """

        (initial, accepting) = self._thompson()
        accepting.accepting = accepting_id
        return NFA(initial)

    def _thompson(self):
        """
        Thompson's construction: convert this AST to an NFA.

        The temporary representation of the NFA is an (initial state, accepting
        state) tuple. The accepting state's ID is None.

        """

        raise NotImplementedError


class SymbolAST(AST):
    """AST leaf node: symbol in the alphabet.

    Attributes:
    symbol -- The symbol for this node.

    """

    def __init__(self, symbol):
        """Create a new symbol AST node.

        Arguments:
        symbol -- The symbol (character) for this node.

        """

        super().__init__()
        self.symbol = symbol

    def _thompson(self):
        initial = NFAState()
        accepting = NFAState()
        initial.add_transition(self.symbol, accepting)
        return (initial, accepting)

    def __repr__(self):
        return 'SymbolAST({})'.format(repr(self.symbol))


class KleeneAST(AST):
    """AST node for the Kleene star operator.

    Attributes:
    operand -- Operand of the closure.

    """

    def __init__(self, operand):
        """Create a new Kleene closure AST node.

        Arguments:
        operand -- AST operand of the closure.

        """

        super().__init__()
        self.operand = operand

    def _thompson(self):
        (initial, accepting) = self.operand._thompson()

        initial.add_transition(None, accepting)
        accepting.add_transition(None, initial)

        return (initial, accepting)

    def __repr__(self):
        return 'KleeneAST({})'.format(repr(self.operand))


class PositiveAST(AST):
    """AST node for the positive closure.

    Attributes:
    operand -- Operand of the closure.

    """

    def __init__(self, operand):
        """Create a new positive closure AST node.

        Arguments:
        operand -- AST operand of the closure.

        """

        super().__init__()
        self.operand = operand

    def _thompson(self):
        (initial, accepting) = self.operand._thompson()

        accepting.add_transition(None, initial)

        return (initial, accepting)

    def __repr__(self):
        return 'PositiveAST({})'.format(repr(self.operand))


class AlternationAST(AST):
    """Alternation (a.k.a. union) of two or more regular expressions.

    Attributes:
    operands -- Tuple containing the alternated regular expressions.

    """

    def __init__(self, *operands):
        """Create a new alternation AST node.

        Arguments:
        operands -- Two or more children AST nodes.

        """

        super().__init__()

        if len(operands) < 2:
            raise ValueError('must alternate two or more operands')

        self.operands = ()
        for ast in operands:
            if isinstance(ast, AlternationAST):
                self.operands += ast.operands
            else:
                self.operands += (ast,)

    def _thompson(self):
        initial = NFAState()
        accepting = NFAState()

        for ast in self.operands:
            (alternate_initial, alternate_accepting) = ast._thompson()
            initial.add_transition(None, alternate_initial)
            alternate_accepting.add_transition(None, accepting)

        return (initial, accepting)

    def __repr__(self):
        return 'AlternationAST({})'.format(', '.join(repr(o) for o in self.operands))


class ConcatenationAST(AST):
    """Concatenation of two or more regular expressions.

    Attributes:
    operands -- Tuple containing the concatenated regular expressions.

    """

    def __init__(self, *operands):
        """Create a new concatenation AST node.

        Arguments:
        operands -- Two or more children AST nodes.

        """

        super().__init__()

        if len(operands) < 2:
            raise ValueError('must concatenate two or more operands')

        self.operands = ()
        for ast in operands:
            if isinstance(ast, ConcatenationAST):
                self.operands += ast.operands
            else:
                self.operands += (ast,)

    def _thompson(self):
        (initial, accepting) = self.operands[0]._thompson()

        for i in range(1, len(self.operands)):
            (next_initial, next_accepting) = self.operands[i]._thompson()
            accepting.add_transition(None, next_initial)
            accepting = next_accepting

        return (initial, accepting)

    def __repr__(self):
        return 'ConcatenationAST({})'.format(', '.join(repr(o) for o in self.operands))


def asts_to_nfa(asts):
    """Convert a list of ASTs to an NFA.

    Thompson's construction is applied to each AST and an initial state is
    created with epsilon transitions to the initial transitions of each
    constructed NFA. The accepting states are given unique IDs ascending from 1
    in the original order of the list.

    """

    initial = NFAState()

    for i, ast in enumerate(asts, 1):
        (ainitial, aaccepting) = ast._thompson()
        aaccepting.accepting = i
        initial.add_transition(None, ainitial)

    return NFA(initial)
