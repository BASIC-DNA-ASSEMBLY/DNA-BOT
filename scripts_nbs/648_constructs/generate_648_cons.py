class BasicCon:
    def __init__(self):
        


class MakeCons:
    @property
    def cons(self):
        return self._cons
    
    @cons.setter
    def cons(self, *basic_cons):
        if not all(isinstance(basic_con, BasicCon) for basic_con in basic_cons):
            raise TypeError("Some basic_cons args not type class BasicCon")
        return self._cons = basic_cons

    
