class Cell:

    #diagonal, top and left
    def __init__(self,score=0,diagonal=0,top=0,left=0):
        self.score=score
        self.diagonal=diagonal
        self.top=top
        self.left=left

    def get_score(self):
        return self.score

    def get_diagonal(self):
        return self.diagonal

    def get_top(self):
        return self.top

    def get_left(self):
        return self.left

    def set_score(self, score):
        self.score = score

    def set_diagonal(self, diagonal):
        self.diagonal = diagonal

    def set_top(self, top):
        self.top = top

    def set_left(self, left):
        self.left = left