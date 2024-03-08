from pysmt.shortcuts import *


class FormatCFA():

    def __init__(self, version):
        self.version = version
    
    @property
    def format_cfa_function(self):
        if self.version in ['_ver2', '_ver3', '_ver3_to_link_B']:
            f = eval( f'self.format_CFA_for_Multiplier{self.version}_in_pysmt')
        elif self.version == '_ver2_to_link_B':
            f = self._format_controlled_Fulladder_in_pysmt
        elif self.version == '_ver1_to_link_B':
            f = self._format_Fulladder_in_pysmt
        else:
            f = None
        return f

    def _format_Fulladder_in_pysmt(self, ws):
        in2, in1, c_in = ws['in2'], ws['in1'], ws['c_in']
        c_out_result = Or(And(c_in, Or(in1, in2)), And(in1, in2))
        out_result = Xor(Xor(in1, in2), c_in)
        return And(Iff(c_out_result, ws['c_out']), Iff(out_result, ws['out']))

    def _format_controlled_Fulladder_in_pysmt(self, vs):
        ws = {v: x for v, x in vs.items() if v in ('in2', 'in1', 'c_in', 'c_out', 'out')}
        ws['in1'] = And(vs['enable'], vs['in1'])
        return self._format_Fulladder_in_pysmt(ws)

    def format_CFA_for_Multiplier_ver2_in_pysmt(self, vs):
        return And(self._format_controlled_Fulladder_in_pysmt(vs), Iff(vs['enable_out'], vs['enable']))

    def format_CFA_for_Multiplier_ver3_in_pysmt(self, vs):
        return And(self.format_CFA_for_Multiplier_ver2_in_pysmt(vs), Iff(vs['in1_out'], vs['in1']))

    def format_CFA_for_Multiplier_ver3_to_link_B_in_pysmt(self, vs):
        return And(self._format_controlled_Fulladder_in_pysmt(), Iff(vs['in1_out'], vs['in1']))

