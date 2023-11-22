import numpy as np
import random as rd
from itertools import combinations

######################################################################################################################################################################
#                                                                                                                                                                    #
#                                                                     Classe                                                                                         #
#                                                                                                                                                                    #
######################################################################################################################################################################


class Particula:
    """Define physics of elastic colisao."""
    
    def __init__(self, massa, raio, posicao, velocidade, tipo = 'atomo'):
        """Initialize a Particle object
        
        massa the massa of particle
        raio the raio of particle
        posicao the posicao vector of particle
        velocidade the velocidade vector of particle
        """
        self.tipo = tipo
        self.massa = massa
        self.raio = raio
        
        # last posicao and velocidade
        self.posicao = np.array(posicao)
        self.velocidade = np.array(velocidade)
        
        # all posicao and velocities recorded during the simulation
        self.todas_posicoes = [np.copy(self.posicao)]
        self.todas_velocidades = [np.copy(self.velocidade)]
        self.todas_magnitudesv = [np.linalg.norm(np.copy(self.velocidade))]
        
    def prox_passo(self, passo):
        """Compute posicao of next passo."""
        self.posicao += passo * self.velocidade
        self.todas_posicoes.append(np.copy(self.posicao)) 
        self.todas_velocidades.append(np.copy(self.velocidade)) 
        self.todas_magnitudesv.append(np.linalg.norm(np.copy(self.velocidade))) 
    
    def checar_colisão(self, particle):
        """Check if there is a colisao with another particle."""
        
        r1, r2 = self.raio, particle.raio
        x1, x2 = self.posicao, particle.posicao
        di = x2-x1
        norm = np.linalg.norm(di)

        if norm-(r1+r2)*1.1 < 0:
            return True
        else:
            return False
        
#         return np.hypot(*(self.raio*2-particle.raio*2)) < self.raio*2 + particle.raio*2


    def realiza_colisao(self, particle, passo):
        """Compute velocidade after colisao with another particle."""
        m1, m2 = self.massa, particle.massa
        r1, r2 = self.raio, particle.raio
        v1, v2 = self.velocidade, particle.velocidade
        x1, x2 = self.posicao, particle.posicao
        di = x1-x2 # mudamos para ficar igual a fórmula
        norm = np.linalg.norm(di)
        if norm-(r1+r2)*1.1 < passo*abs(np.dot(v1-v2, di))/norm:
            self.velocidade = v1 - 2. * m2/(m1+m2) * np.dot(v1-v2, di) / (np.linalg.norm(di)**2.) * di
            particle.velocidade = v2 - 2. * m1/(m2+m1) * np.dot(v2-v1, (-di)) / (np.linalg.norm(di)**2.) * (-di)
            

    def colisao_paredess(self, passo, size, probabilidade_reacao_paredes):
        """Compute velocidade after hitting an edge.
        passo the computation passo
        size the medium size
        """
        r, v, x = self.raio, self.velocidade, self.posicao
        projx = passo*abs(np.dot(v,np.array([1.,0.])))
        projy = passo*abs(np.dot(v,np.array([0.,1.])))
        if abs(x[0])-r < projx or abs(size-x[0])-r < projx:
            self.velocidade[0] *= -1
            probabilidade_reacao_paredes = 0.98

            
        if abs(x[1])-r < projy or abs(size-x[1])-r < projy:
            self.velocidade[1] *= -1.
            probabilidade_reacao_paredes = 0.98  

        return probabilidade_reacao_paredes

            
       
            
######################################################################################################################################################
#                                                                                                                                                    #
#                                                                    Funções                                                                         #
#                                                                                                                                                    #
######################################################################################################################################################

def mudar_passo(lista_particulas, passo, size, probabilidade_reacao):
    """Solve a passo for every particle."""
    
    # Detect edge-hitting and colisao of every particle
    for i in range(len(lista_particulas)):
        lista_particulas[i].colisao_paredess(passo,size, probabilidade_reacao)
        for j in range(i+1,len(lista_particulas)):
                lista_particulas[i].realiza_colisao(lista_particulas[j],passo)    

                
    # Computa a posicao de todas particulas  
    for particula in lista_particulas:
        particula.prox_passo(passo)
        
        probabilidade_reacao = particula.colisao_paredess(passo, size, probabilidade_reacao)




def gerar_particulas(N, raio, massa, tamanho_caixa, tipo):
    """Generate N Particle objects in a random way in a list."""
    lista_particulas = []

    for i in range(N):
        
        magnitude_velocidade = np.random.rand(1)*25
        angulo_velocidade = np.random.rand(1)*2*np.pi
        v = np.append(magnitude_velocidade*np.cos(angulo_velocidade), magnitude_velocidade*np.sin(angulo_velocidade))
        
        colisao = True
        while(colisao == True):
            
            colisao = False
            pos = raio + np.random.rand(2)*(tamanho_caixa-2*raio) 
            nova_particula = Particula(massa, raio, pos, v, tipo)
            for j in range(len(lista_particulas)):

                colisao = nova_particula.checar_colisão( lista_particulas[j] )

                if colisao == True:
                    break

        lista_particulas.append(nova_particula)
    return lista_particulas

#########################################################################################################################################################
#                                                                                                                                                       #
#                                                                    Reação                                                                             #
#                                                                                                                                                       #
#########################################################################################################################################################


def simular_reacao(lista_particulas, probabilidade_reacao):    
    for particula1, particula2 in combinations(lista_particulas, 2):
                
        if particula1.checar_colisão(particula2) and particula1.tipo == 'atomo' and particula2.tipo == 'atomo':

            valor_aleatorio = 0

            if valor_aleatorio < probabilidade_reacao:
                nova_massa = particula1.massa + particula2.massa
                nova_posicao = (particula1.posicao + particula2.posicao) / 2 #centro de massa
                nova_velocidade = (particula1.massa * particula1.velocidade + particula2.massa * particula2.velocidade) / nova_massa
                nova_particula = Particula(
                    massa=nova_massa,
                    raio=(particula1.raio**3 + particula2.raio**3)**(1/3),  # Lei da conservação de volume
                    tipo='molecula',
                    posicao=nova_posicao,
                    velocidade=nova_velocidade
                ) 
                
                lista_particulas.remove(particula1)
                lista_particulas.remove(particula2)
                lista_particulas.append(nova_particula)
                return nova_particula
            else:
                #colisão elastica
                return realiza_colisao(particula1, particula2, passo)
    else:
        # Não houve colisão
        return None   

#########################################################################################################################################################
#                                                                                                                                                       #
#                                                                    Desafio 3                                                                          #
#                                                                                                                                                       #
#########################################################################################################################################################



def gerar_particulas_dois_sistemas(N, raio, massa, tamanho_caixa, tipo):
    """Generate N Particle objects in a random way in a list."""
    lista_particulas_sistema1, lista_particulas_sistema2, lista_particulas_sistema3, lista_particulas_sistema4, lista_particulas_sistema5 = [], [], [], [], []

    for i in range(N):
        
        magnitude_velocidade_sistema1 = np.random.rand(1)*25  # Velocidade para o primeiro sistema
        magnitude_velocidade_sistema2 = np.random.rand(1)*50  # Velocidade para o segundo sistema
        magnitude_velocidade_sistema3 = np.random.rand(1)*75  # Velocidade para o terceiro sistema
        magnitude_velocidade_sistema4 = np.random.rand(1)*100  # Velocidade para o quarto sistema
        magnitude_velocidade_sistema5 = np.random.rand(1)*125  # Velocidade para o quinto sistema
        
        
        angulo_velocidade = np.random.rand(1)*2*np.pi
        
        v_sistema1 = np.append(magnitude_velocidade_sistema1 * np.cos(angulo_velocidade), magnitude_velocidade_sistema1 * np.sin(angulo_velocidade))
        v_sistema2 = np.append(magnitude_velocidade_sistema2 * np.cos(angulo_velocidade), magnitude_velocidade_sistema2 * np.sin(angulo_velocidade))
        v_sistema3 = np.append(magnitude_velocidade_sistema3 * np.cos(angulo_velocidade), magnitude_velocidade_sistema3 * np.sin(angulo_velocidade))
        v_sistema4 = np.append(magnitude_velocidade_sistema4 * np.cos(angulo_velocidade), magnitude_velocidade_sistema4 * np.sin(angulo_velocidade))
        v_sistema5 = np.append(magnitude_velocidade_sistema5 * np.cos(angulo_velocidade), magnitude_velocidade_sistema5 * np.sin(angulo_velocidade))

        
        pos = raio + np.random.rand(2)*(tamanho_caixa-2*raio) 
        
        nova_particula_sistema1 = Particula(massa, raio, pos, v_sistema1, tipo)
        nova_particula_sistema2 = Particula(massa, raio, pos, v_sistema2, tipo)
        nova_particula_sistema3 = Particula(massa, raio, pos, v_sistema3, tipo)
        nova_particula_sistema4 = Particula(massa, raio, pos, v_sistema4, tipo)
        nova_particula_sistema5 = Particula(massa, raio, pos, v_sistema5, tipo)

        colisao_sistema1, colisao_sistema2, colisao_sistema3, colisao_sistema4, colisao_sistema5 = False, False, False, False, False
        
        for particula_sistema1 in lista_particulas_sistema1:
            colisao_sistema1 = nova_particula_sistema1.checar_colisão(particula_sistema1)
            if colisao_sistema1:
                break
                
        for particula_sistema2 in lista_particulas_sistema2:
            colisao_sistema2 = nova_particula_sistema2.checar_colisão(particula_sistema2)
            if colisao_sistema2:
                break
                
        for particula_sistema3 in lista_particulas_sistema3:
            colisao_sistema3 = nova_particula_sistema3.checar_colisão(particula_sistema3)
            if colisao_sistema3:
                break

        for particula_sistema4 in lista_particulas_sistema4:
            colisao_sistema4 = nova_particula_sistema4.checar_colisão(particula_sistema4)
            if colisao_sistema4:
                break
                
        for particula_sistema5 in lista_particulas_sistema5:
            colisao_sistema5 = nova_particula_sistema5.checar_colisão(particula_sistema5)
            if colisao_sistema5:
                break

        if not colisao_sistema1:
            lista_particulas_sistema1.append(nova_particula_sistema1)
            
        if not colisao_sistema2:
            lista_particulas_sistema2.append(nova_particula_sistema2)
            
        if not colisao_sistema3:
            lista_particulas_sistema3.append(nova_particula_sistema3)
        
        if not colisao_sistema4:
            lista_particulas_sistema4.append(nova_particula_sistema4)
            
        if not colisao_sistema5:
            lista_particulas_sistema5.append(nova_particula_sistema5)
        
    return lista_particulas_sistema1, lista_particulas_sistema2, lista_particulas_sistema3, lista_particulas_sistema4, lista_particulas_sistema5



def exponencial(t, a, k, c):
    return a * np.exp(-k * t) + c



#########################################################################################################################################################
#                                                                                                                                                       #
#                                                                    Catalisador                                                                        #
#                                                                                                                                                       #
#########################################################################################################################################################




# Classe Particula_Catalisador
class Particula_Catalisador:
    def __init__(self, massa, raio, posicao, velocidade, tipo='parede'):
        self.tipo = tipo
        self.massa = massa
        self.raio = raio
        self.posicao = np.array(posicao)
        self.velocidade = np.array(velocidade)

    def checar_colisao(self, particle):
        r1, r2 = self.raio, particle.raio
        x1, x2 = self.posicao, particle.posicao
        di = x2 - x1
        norm = np.linalg.norm(di)

        if norm - (r1 + r2) * 1.1 < 0:
            return True
        else:
            return False

    def realiza_colisao(self, particle, passo):
        m1, m2 = self.massa, particle.massa
        r1, r2 = self.raio, particle.raio
        v1, v2 = self.velocidade, particle.velocidade
        x1, x2 = self.posicao, particle.posicao
        di = x1 - x2
        norm = np.linalg.norm(di)
        
        if norm - (r1 + r2) * 1.1 < passo * abs(np.dot(v1 - v2, di)) / norm:
            self.velocidade = v1 - 2. * m2 / (m1 + m2) * np.dot(v1 - v2, di) / (np.linalg.norm(di) ** 2.) * di
            particle.velocidade = v2 - 2. * m1 / (m2 + m1) * np.dot(v2 - v1, (-di)) / (np.linalg.norm(di) ** 2.) * (-di)

    def colisao_paredes(self, passo, size):
        r, v, x = self.raio, self.velocidade, self.posicao
        projx = passo * abs(np.dot(v, np.array([1., 0.])))
        projy = passo * abs(np.dot(v, np.array([0., 1.])))
        if abs(x[0]) - r < projx or abs(size - x[0]) - r < projx:
            self.velocidade[0] *= -1
        if abs(x[1]) - r < projy or abs(size - x[1]) - r < projy:
            self.velocidade[1] *= -1.

# Função para verificar colisões com a parede
def verifica_colisao_parede(particula, size):
    r = particula.raio
    x = particula.posicao
    v = particula.velocidade
    projx = abs(np.dot(v, np.array([1., 0.])))
    projy = abs(np.dot(v, np.array([0., 1.])))

    if abs(x[0]) - r < projx or abs(size - x[0]) - r < projx:
        return True
    if abs(x[1]) - r < projy or abs(size - x[1]) - r < projy:
        return True
    return False

# Função para modificar a probabilidade de reação se a partícula colidir com a parede
def mudar_probabilidade(particula, probabilidade, taxa_aumento):
    if verifica_colisao_parede(particula, size):
        probabilidade *= taxa_aumento
    return probabilidade  